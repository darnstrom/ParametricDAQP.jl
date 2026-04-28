using Printf 

function write_array(f,A,name,type)
    N = length(A)
    write(f,"$type $name[$N] = {\n")
    for i in 1:N
        write(f, "($type)$(A[i]),\n")
    end
    write(f,"};\n")
end

function codegen(sol::Solution;dir="codegen",fname="pdaqp", float_type="float", c_float_store=float_type, int_type="unsigned short",
        max_reals=1e12, dual = false, bfs=true, clipping=false,
        store_transpose=false, store_offset=false)
    bst = build_tree(sol;max_reals,dual,bfs,clipping);
    isnothing(bst) && return -1
    codegen(bst;dir,fname,float_type,int_type,c_float_store,store_transpose,store_offset)
    return 1
end

function codegen(bst::BinarySearchTree; dir="codegen",fname="pdaqp", float_type="float", c_float_store=float_type, int_type="unsigned short", store_transpose=false, store_offset=false)
    isdir(dir) || mkdir(dir)
    # Get number of outputs
    nth,nz = size(bst.feedbacks[1]).-(1,0)
    # Concatenate feedbacks into one array (normal layout)
    feedbacks = reduce(hcat,bst.feedbacks)
    # Optionally also build transposed layout
    feedbacks_T = store_transpose ? reduce(hcat,[collect(f') for f in bst.feedbacks]) : nothing

    if(!isempty(bst.duals))
        duals = reduce(hcat,bst.duals)
    else
        duals =  zeros(0,0)
    end

    # Write header file
    fh = open(joinpath(dir,fname*".h"), "w")
    hguard = uppercase(fname)*"_H"
    write(fh, "#ifndef $hguard\n")
    write(fh, "#define $hguard\n\n")

    write(fh, "typedef $float_type c_float;\n")
    write(fh, "typedef $c_float_store c_float_store;\n")
    write(fh, "typedef $int_type c_int;\n")
    write(fh, "#define $(uppercase(fname))_N_PARAMETER $nth\n")
    write(fh, "#define $(uppercase(fname))_N_SOLUTION $nz\n\n")
    !isempty(duals) && write(fh, "#define $(uppercase(fname))_N_CONSTRAINTS $(size(bst.duals[1],2))\n\n")
    eval_sol_args = isempty(duals) ? "c_float* solution" : "c_float* solution, c_float* dual"
    write(fh, "void $(fname)_evaluate(c_float* parameter, $eval_sol_args);\n")
    write(fh, "extern c_float_store $(fname)_feedbacks[$(length(feedbacks))];\n")
    store_transpose && write(fh, "extern c_float_store $(fname)_feedbacks_T[$(length(feedbacks))];\n")
    store_offset && write(fh, "extern volatile int $(fname)_feedback_offset;\n")

    !isempty(bst.clipping) && write(fh, "c_float $(fname)_clip(c_float v, c_float min, c_float max);\n")
    write(fh, "#endif // ifndef $hguard\n");
    close(fh)

    # Write source file
    fsrc = open(joinpath(dir,fname*".c"), "w")
    write(fsrc, "#include \"$fname.h\"\n")
    write_array(fsrc,bst.halfplanes,fname*"_halfplanes","c_float_store")
    write_array(fsrc,feedbacks,fname*"_feedbacks","c_float_store")
    store_offset && write(fsrc, "volatile int $(fname)_feedback_offset;\n");
    !isnothing(feedbacks_T) && write_array(fsrc,feedbacks_T,fname*"_feedbacks_T","c_float_store")
    !isempty(duals) && write_array(fsrc,duals,fname*"_duals","c_float")
    !isempty(bst.clipping) && write_array(fsrc,bst.clipping[:,1],fname*"_out_min","c_float")
    !isempty(bst.clipping) && write_array(fsrc,bst.clipping[:,2],fname*"_out_max","c_float")
    # -1 to indices since 0 index C is 1 index Julia
    write_array(fsrc,bst.hp_list.-1,fname*"_hp_list","c_int")
    write_array(fsrc,bst.jump_list,fname*"_jump_list","c_int")

    clip_call = isempty(bst.clipping) ? "val" : "$(fname)_clip(val,$(fname)_out_min[i],$(fname)_out_max[i])"
    clip_call_t = isempty(bst.clipping) ? "solution[i]" : "$(fname)_clip(solution[i],$(fname)_out_min[i],$(fname)_out_max[i])"
    src_code = """void $(fname)_evaluate(c_float* parameter, $eval_sol_args){
        int i,j,disp;
        int id,next_id;
        c_float val;
        id = 0;
        next_id = id+$(fname)_jump_list[id];
        while(next_id != id){
            // Compute halfplane value
            disp = $(fname)_hp_list[id]*($(uppercase(fname))_N_PARAMETER+1);
            for(i=0, val=0; i<$(uppercase(fname))_N_PARAMETER; i++)
                val += parameter[i] * $(fname)_halfplanes[disp++];
            id = next_id + (val <= $(fname)_halfplanes[disp]);
            next_id = id+$(fname)_jump_list[id];
        }
        // Leaf node reached -> evaluate affine function
        disp = $(fname)_hp_list[id]*($(uppercase(fname))_N_PARAMETER+1)*$(uppercase(fname))_N_SOLUTION;
        for(i=0; i < $(uppercase(fname))_N_SOLUTION; i++){
            for(j=0, val=0; j < $(uppercase(fname))_N_PARAMETER; j++)
                val += parameter[j] * $(fname)_feedbacks[disp++];
            val += $(fname)_feedbacks[disp++];
            solution[i] = $clip_call;
        }
    """
    if store_offset
        src_code *= """    $(fname)_feedback_offset = $(fname)_hp_list[id]*($(uppercase(fname))_N_PARAMETER+1)*$(uppercase(fname))_N_SOLUTION;
        """
    end
    if(isempty(duals))
        src_code *= """
        }
        """
    else
        src_code *= """    disp = $(fname)_hp_list[id]*($(uppercase(fname))_N_PARAMETER+1)*$(uppercase(fname))_N_CONSTRAINTS;
            for(i=0; i < $(uppercase(fname))_N_CONSTRAINTS; i++){
                for(j=0, val=0; j < $(uppercase(fname))_N_PARAMETER; j++)
                    val += parameter[j] * $(fname)_duals[disp++];
                val += $(fname)_duals[disp++];
                dual[i] = val;
            }
        }
        """
    end

    if(!isempty(bst.clipping))
        src_code *= """
        c_float $(fname)_clip(c_float v, c_float min, c_float max) {
            const c_float vmin = v < min ? min : v;
            return vmin > max ? max : vmin;
        }
        """
    end

    write(fsrc,src_code);
    close(fsrc)

    # Write simple example
    fex = open(joinpath(dir,"example.c"), "w")
    dual_array = isempty(duals) ? "" : "c_float dual[$(size(bst.duals[1],2))];"
    dual_arg = isempty(duals) ? "" : ",dual"
    write(fex, """
#include "$(fname).h"
#include <stdio.h>

int main(){
    c_float solution[$nz];
    c_float parameter[$nth];
    $dual_array

    int i;
    // Initialize parameter
    for(i=0; i< $nth; i++)
        parameter[i] = 0;

    // Get the solution at the parameter
    $(fname)_evaluate(parameter,solution$dual_arg);

    printf("For the parameter\\n");
    for(i=0; i< $nth; i++)
        printf("%f\\n",parameter[i]);
    printf("the solution is\\n");
    for(i=0; i< $nz; i++)
        printf("%f\\n",solution[i]);
}
          """)
    close(fex)
end

function codegen_implicit(mpp;fname="pdaqp_workspace", dir="codegen", opt_settings=nothing, src=true, float_type="double",warm_start=false)
    length(dir)==0 && (dir="codegen")
    dir[end] != '/' && (dir*="/") ## Make sure it is a correct directory path

    output_transform = get_codegen_output_transform(mpp)

    # Generate mpldp
    mpldp = setup_mpp(mpp)
    mpldp isa MPLDP || error("codegen_daqp requires a positive definite Hessian H")
    size(mpldp.d, 1) == mpldp.n_theta + 2 || error("codegen_implicit requires a problem with upper and lower affine bounds")

    # Generate DAQP workspace
    d = DAQPBase.Model() 
    !isnothing(opt_settings) && DAQPBase.settings(d,opt_settings)
    m = length(mpp.bu)
    senses  = hasproperty(mpp, :senses) ? Vector{Cint}(mpp.senses) : zeros(Cint, m)
    DAQPBase.setup(d, Matrix{Cdouble}(mpp.H), vec(Cdouble.(mpp.f)),
                   Matrix{Cdouble}(mpp.A), vec(Cdouble.(mpp.bu)),
                   vec(Cdouble.(mpp.bl)), senses)

    DAQPBase.codegen(d;fname,dir,src=false)
    src && copy_daqp_codegen_headers(dir)

    if( float_type == "float" || float_type == "single")
        # Append #define DAQP_SINGLE_PRECISION at the top of types
        mv(joinpath(dir,"types.h"),joinpath(dir,"types_old.h"))
        fold = open(joinpath(dir,"types_old.h"),"r")
        s = read(fold, String)
        fnew = open(joinpath(dir,"types.h"),"w")
        write(fnew, "#ifndef DAQP_SINGLE_PRECISION\n # define DAQP_SINGLE_PRECISION\n#endif \n"*s);
        close(fold)
        close(fnew)
        rm(joinpath(dir,"types_old.h"))
    end

    # Append MPP-specific data/functions
    render_pdaqp_workspace(mpldp;fname,dir,float_type, fmode="a",warm_start, output_transform)
    @info "Generated code for parameteric program" dir fname
end

function copy_daqp_codegen_headers(dir)
    include_dir = normpath(dirname(DAQPBase.libdaqp), "..", "include", "daqp")
    for file in ("daqp.h", "auxiliary.h", "factorization.h", "constants.h", "types.h")
        cp(joinpath(include_dir, file), joinpath(dir, file); force=true)
    end
end

function get_codegen_output_transform(mpp)
    f_theta = hasproperty(mpp, :f_theta) ? mpp.f_theta : mpp.F
    nth = size(f_theta, 2)
    n = size(mpp.H, 1)

    if hasproperty(mpp, :post_transform) && !isnothing(mpp.post_transform)
        Z, affine = mpp.post_transform
        Z = Matrix{Float64}(Z)
        affine = Matrix{Float64}(affine)
        size(Z, 1) == size(affine, 1) || error("post_transform components must have the same number of rows")
        size(affine, 2) == nth + 1 || error("post_transform affine term must have nth+1 columns")
        out_inds = hasproperty(mpp, :out_inds) && !isnothing(mpp.out_inds) ? collect(mpp.out_inds) : collect(1:size(Z, 1))
    else
        Z = Matrix{Float64}(I, n, n)
        affine = zeros(Float64, n, nth + 1)
        out_inds = hasproperty(mpp, :out_inds) && !isnothing(mpp.out_inds) ? collect(mpp.out_inds) : collect(1:n)
    end

    maximum(out_inds) <= size(Z, 1) || error("out_inds contains an index outside the post-transform output")
    return Z[out_inds, :], affine[out_inds, :]
end

function render_pdaqp_workspace(mpldp;fname="pdaqp_workspace",dir="",fmode="w", float_type="double", warm_start=false, output_transform)
    nth = mpldp.n_theta
    n_reduced = size(output_transform[1], 2)
    n_output = size(output_transform[1], 1)
    m = length(mpldp.norm_factors)
    Dth = mpldp.d[1:nth, :]
    du = mpldp.d[nth + 1, :]
    dl = mpldp.d[nth + 2, :]
    Z_map, affine_offset = output_transform

    # Setup files
    fh = open(dir*fname*".h", fmode)
    fsrc = open(dir*fname*".c", fmode)

    # HEADER 
    hguard = uppercase(fname)*"_PDAQP_H"
    @printf(fh, "#ifndef %s\n",   hguard);
    @printf(fh, "#define %s\n\n", hguard);

    @printf(fh, "#define PDAQP_N_PARAMETERS %d\n",nth);
    @printf(fh, "#define PDAQP_N_CONSTRAINTS %d\n",m);
    @printf(fh, "#define PDAQP_N_REDUCED_SOLUTION %d\n",n_reduced);
    @printf(fh, "#define PDAQP_N_SOLUTION %d\n",n_output);

    @printf(fh, "#define DAQP_WARMSTART %d\n\n", warm_start ? 1 : 0)

    @printf(fh, "extern c_float mpqp_parameter[%d];\n", nth);

    @printf(fh, "extern c_float Dth[%d];\n", nth*m);
    @printf(fh, "extern c_float du[%d];\n", m);
    @printf(fh, "extern c_float dl[%d];\n\n", m);

    @printf(fh, "extern c_float solution_map[%d];\n", n_output*n_reduced);
    @printf(fh, "extern c_float solution_offset[%d];\n\n", n_output*(nth+1));
    @printf(fh, "void %s_form_qp(const c_float* parameter);\n", fname);
    @printf(fh, "int %s_solve(const c_float* parameter, c_float* solution);\n\n", fname);


    # SRC 
    write_array(fsrc,zeros(nth),"mpqp_parameter","c_float");
    write_array(fsrc,Dth,"Dth", "c_float");
    write_array(fsrc,du,"du", "c_float");
    write_array(fsrc,dl,"dl", "c_float");
    write_array(fsrc,copy(Z_map'),"solution_map", "c_float");
    write_array(fsrc,copy(affine_offset'),"solution_offset", "c_float");

    @printf(fsrc, "#include \"daqp.h\"\n");
    write(fsrc, """
void $(fname)_form_qp(const c_float* parameter){
    int i,j,disp;
    c_float val;

    for(i = 0; i < $(nth); ++i){
        mpqp_parameter[i] = parameter[i];
    }

    disp = 0;
    for(i = 0; i < $(m); ++i){
        for(j = 0, val = 0; j < $(nth); ++j){
            val += parameter[j] * Dth[disp++];
        }
        daqp_work.dupper[i] = du[i] + val;
        daqp_work.dlower[i] = dl[i] + val;
    }
}

int $(fname)_solve(const c_float* parameter, c_float* solution){
    int i,j,disp,exitflag;
    c_float val;

    $(fname)_form_qp(parameter);
#if !DAQP_WARMSTART
    reset_daqp_workspace(&daqp_work);
#endif
    exitflag = daqp_ldp(&daqp_work);
    if(exitflag > 0){
        ldp2qp_solution(&daqp_work);
        disp = 0;
        for(i = 0; i < $(n_output); ++i){
            for(j = 0, val = 0; j < $(n_reduced); ++j){
                val += solution_map[disp++] * daqp_work.x[j];
            }
            for(j = 0; j < $(nth); ++j){
                val += solution_offset[i * ($(nth) + 1) + j] * parameter[j];
            }
            val += solution_offset[(i + 1) * ($(nth) + 1) - 1];
            solution[i] = val;
        }
    }
    return exitflag;
}
""")

    @printf(fh, "#endif // ifndef %s\n", hguard);

    close(fh)
    close(fsrc)
end
