function write_array(f,A,name,type)
    N = length(A)
    write(f,"$type $name[$N] = {\n")
    for i in 1:N 
        write(f, "($type)$(A[i]),\n")
    end
    write(f,"};\n")
end

function codegen(sol::Solution;dir="codegen",fname="pdaqp", float_type="float", int_type="unsigned short",
        max_reals=1e12, dual = false)
    bst = build_tree(sol;max_reals,dual);
    isnothing(bst) && return -1
    codegen(bst;dir,fname,float_type,int_type)
    return 1
end

function codegen(bst::BinarySearchTree; dir="codegen",fname="pdaqp", float_type="float", int_type="unsigned short")
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
    write(fh, "typedef $int_type c_int;\n")
    write(fh, "#define $(uppercase(fname))_N_PARAMETER $nth\n")
    write(fh, "#define $(uppercase(fname))_N_SOLUTION $nz\n\n")
    !isempty(duals) && write(fh, "#define $(uppercase(fname))_N_CONSTRAINTS $(size(bst.duals[1],2))\n\n")
    eval_sol_args = isempty(duals) ? "c_float* solution" : "c_float* solution, c_float* dual"
    write(fh, "void $(fname)_evaluate(c_float* parameter, $eval_sol_args);\n")
    write(fh, "#endif // ifndef $hguard\n");
    close(fh)

    # Write source file
    fsrc = open(joinpath(dir,fname*".c"), "w")
    write(fsrc, "#include \"$fname.h\"\n")
    write_array(fsrc,bst.halfplanes,fname*"_halfplanes","c_float")
    write_array(fsrc,feedbacks,fname*"_feedbacks","c_float")
    !isempty(duals) && write_array(fsrc,duals,fname*"_duals","c_float")
    # -1 to indices since 0 index C is 1 index Julia
    write_array(fsrc,bst.hp_list.-1,fname*"_hp_list","c_int")
    write_array(fsrc,bst.jump_list,fname*"_jump_list","c_int") 

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
            solution[i] = val;
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
