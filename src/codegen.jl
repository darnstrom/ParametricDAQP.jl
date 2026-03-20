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
        max_reals=1e12, dual = false, bfs=true, clipping=false)
    bst = build_tree(sol;max_reals,dual,bfs,clipping);
    isnothing(bst) && return -1
    codegen(bst;dir,fname,float_type,int_type,c_float_store)
    return 1
end

function codegen(bst::BinarySearchTree; dir="codegen",fname="pdaqp", float_type="float", c_float_store=float_type, int_type="unsigned short")
    isdir(dir) || mkdir(dir)
    # Get number of outputs 
    nth,nz = size(bst.feedbacks[1]).-(1,0)
    # Concatenate feedbacks into one array
    feedbacks = reduce(hcat,bst.feedbacks)

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
    !isempty(bst.clipping) && write(fh, "c_float $(fname)_clip(c_float v, c_float min, c_float max);\n")
    write(fh, "#endif // ifndef $hguard\n");
    close(fh)

    # Write source file
    fsrc = open(joinpath(dir,fname*".c"), "w")
    write(fsrc, "#include \"$fname.h\"\n")
    write_array(fsrc,bst.halfplanes,fname*"_halfplanes","c_float_store")
    write_array(fsrc,feedbacks,fname*"_feedbacks","c_float_store")
    !isempty(duals) && write_array(fsrc,duals,fname*"_duals","c_float")
    !isempty(bst.clipping) && write_array(fsrc,bst.clipping[:,1],fname*"_out_min","c_float")
    !isempty(bst.clipping) && write_array(fsrc,bst.clipping[:,2],fname*"_out_max","c_float")
    # -1 to indices since 0 index C is 1 index Julia
    write_array(fsrc,bst.hp_list.-1,fname*"_hp_list","c_int")
    write_array(fsrc,bst.jump_list,fname*"_jump_list","c_int") 

    clip_call = isempty(bst.clipping) ? "val" : "$(fname)_clip(val,$(fname)_out_min[i],$(fname)_out_max[i])"
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

    # Generate mpldp
    mpldp = setup_mpp(mpp)
    mpldp isa MPLDP || error("codegen_daqp requires a positive definite Hessian H")

    # Generate DAQP workspace
    d = DAQPBase.Model() 
    !isnothing(opt_settings) && DAQPBase.settings(d,opt_settings)
    m = size(mpp.bu, 1)
    blower  = fill(-1e30, m)
    senses  = hasproperty(mpp, :senses) ? Vector{Cint}(mpp.senses) : zeros(Cint, m)
    DAQPBase.setup(d, Matrix{Cdouble}(mpp.H), Vector{Cdouble}(mpp.f),
                   Matrix{Cdouble}(mpp.A), Vector{Cdouble}(mpp.bu),
                   Vector{Cdouble}(mpp.bl), senses)

    DAQPBase.codegen(d;fname,dir,src)

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
    render_pdaqp_workspace(mpldp;fname,dir,float_type, fmode="a",warm_start)
    @info "Generated code for parameteric program" dir fname
end

function render_pdaqp_workspace(mpldp;fname="pdaqp_workspace",dir="",fmode="w", float_type="double", warm_start=false)
    nth = mpldp.n_theta
    n_out = size(mpldp.HinvF, 2)  # number of output (decision) variables
    m = length(mpldp.norm_factors)

    # Setup files
    fh = open(dir*fname*".h", fmode)
    fsrc = open(dir*fname*".c", fmode)

    # HEADER
    hguard = uppercase(fname)*"_PDAQP_H"
    @printf(fh, "#ifndef %s\n",   hguard);
    @printf(fh, "#define %s\n\n", hguard);

    @printf(fh, "#define PDAQP_N_PARAMETERS %d\n",  nth);
    @printf(fh, "#define PDAQP_N_CONSTRAINTS %d\n", m);
    @printf(fh, "#define PDAQP_N_DECISION %d\n\n",  n_out);

    if warm_start
        @printf(fh, "#define DAQP_WARMSTART\n\n");
    end

    @printf(fh, "extern c_float mpqp_parameter[%d];\n", nth);
    @printf(fh, "extern c_float Dth[%d];\n", nth*m);
    @printf(fh, "extern c_float du[%d];\n", m);
    @printf(fh, "extern c_float dl[%d];\n\n", m);
    @printf(fh, "extern c_float Z_offset[%d];\n", n_out*nth);
    @printf(fh, "extern c_float z_offset[%d];\n\n", n_out);

    # Function declarations
    @printf(fh, "void pdaqp_update_qp(c_float* th, c_float* dupper, c_float* dlower);\n");
    @printf(fh, "void pdaqp_get_solution(c_float* th, c_float* solution, c_float* xstar);\n");
    @printf(fh, "int  pdaqp_compute_solution(c_float* th, c_float* solution);\n\n");

    @printf(fh, "#endif // ifndef %s\n", hguard);

    close(fh)

    # SRC: global parametric data arrays
    write_array(fsrc, zeros(nth),                       "mpqp_parameter", "c_float");
    write_array(fsrc, mpldp.d[1:end-2,:],               "Dth",            "c_float");
    write_array(fsrc, mpldp.d[end-1,:],                 "du",             "c_float");
    write_array(fsrc, mpldp.d[end,:],                   "dl",             "c_float");
    write_array(fsrc, copy(mpldp.HinvF[1:end-1,:]'),    "Z_offset",       "c_float");
    write_array(fsrc, mpldp.HinvF[end,:],               "z_offset",       "c_float");

    # Include headers (after the data definitions so that include guards prevent
    # double-inclusion of any type definitions already present in the file)
    @printf(fsrc, "#include \"%s.h\"\n", fname);
    @printf(fsrc, "#include \"daqp.h\"\n\n");

    # C function: update DAQP upper/lower bounds given parameter th
    src_code = """void pdaqp_update_qp(c_float* th, c_float* dupper, c_float* dlower){
    int i, j, disp;
    c_float b_shift_th;
    for(i = 0, disp = 0; i < PDAQP_N_CONSTRAINTS; i++){
        b_shift_th = 0;
        for(j = 0; j < PDAQP_N_PARAMETERS; j++) b_shift_th += Dth[disp++]*th[j];
        dupper[i] = du[i] + b_shift_th;
        dlower[i] = dl[i] + b_shift_th;
    }
}

"""
    # C function: transform DAQP QP solution xstar back to the nominal variable
    src_code *= """void pdaqp_get_solution(c_float* th, c_float* solution, c_float* xstar){
    int i, j;
    c_float sol_shift;
    for(i = 0; i < PDAQP_N_DECISION; i++){
        sol_shift = z_offset[i];
        for(j = 0; j < PDAQP_N_PARAMETERS; j++)
            sol_shift += Z_offset[i + j*PDAQP_N_DECISION]*th[j];
        solution[i] = xstar[i] + sol_shift;
    }
}

"""
    # C function: update constraints, solve the LDP, and recover the nominal solution
    src_code *= """int pdaqp_compute_solution(c_float* th, c_float* solution){
    pdaqp_update_qp(th, daqp_work.dupper, daqp_work.dlower);
    daqp_work.reuse_ind = 0;
"""
    if !warm_start
        src_code *= """    daqp_deactivate_constraints(&daqp_work);
    reset_daqp_workspace(&daqp_work);
"""
    end
    src_code *= """    int exitflag = daqp_ldp(&daqp_work);
    ldp2qp_solution(&daqp_work);
    pdaqp_get_solution(th, solution, daqp_work.x);
    return exitflag;
}
"""
    write(fsrc, src_code)
    close(fsrc)
end
