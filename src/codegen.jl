function write_array(f,A,name,type)
    N = length(A)
    write(f,"$type $name[$N] = {\n")
    for i in 1:N 
        write(f, "($type)$(A[i]),\n")
    end
    write(f,"};\n")
end

function codegen(sol::Solution;dir="codegen",fname="pdaqp", float_type="float", int_type="unsigned short")
    bst = build_tree(sol);
    codegen(bst;dir,fname,float_type,int_type)
end

function codegen(bst::BinarySearchTree; dir="codegen",fname="pdaqp", float_type="float", int_type="unsigned short")
    isdir(dir) || mkdir(dir)
    # Get number of outputs 
    nth,nz = size(bst.feedbacks[1]).-(1,0)
    # Concatenate feedbacks into one array
    feedbacks = reduce(hcat,bst.feedbacks)

    # Write header file
    fh = open(joinpath(dir,fname*".h"), "w")
    hguard = uppercase(fname)*"_H"
    write(fh, "#ifndef $hguard\n")
    write(fh, "#define $hguard\n\n")

    write(fh, "typedef $float_type c_float;\n")
    write(fh, "typedef $int_type c_int;\n")
    write(fh, "#define $(uppercase(fname))_N_PARAMETER $nth\n")
    write(fh, "#define $(uppercase(fname))_N_SOLUTION $nz\n\n")
    write(fh, "void $(fname)_evaluate(c_float* parameter, c_float* solution);\n")
    write(fh, "#endif // ifndef $hguard\n");
    close(fh)

    # Write source file
    fsrc = open(joinpath(dir,fname*".c"), "w")
    write(fsrc, "#include \"$fname.h\"\n")
    write_array(fsrc,bst.halfplanes,fname*"_halfplanes","c_float")
    write_array(fsrc,feedbacks,fname*"_feedbacks","c_float")
    # -1 to indices since 0 index C is 1 index Julia
    write_array(fsrc,bst.hp_list.-1,fname*"_hp_list","c_int")
    write_array(fsrc,bst.jump_list.-1,fname*"_jump_list","c_int") 

    write(fsrc, """
void $(fname)_evaluate(c_float* parameter, c_float* solution){
    int i,j,disp;
    int id,next_id;
    c_float val;
    id = 0;
    next_id = $(fname)_jump_list[id];
    while(next_id != 0){
        // Compute halfplane value 
        disp = $(fname)_hp_list[id]*($(uppercase(fname))_N_PARAMETER+1);
        for(i=0, val=0; i<$(uppercase(fname))_N_PARAMETER; i++)
            val += parameter[i] * $(fname)_halfplanes[disp++];
        if(val <= $(fname)_halfplanes[disp])// positive branch
            id = next_id+1;
        else // negative branch
            id = next_id;
        next_id = $(fname)_jump_list[id];
    }
    // Leaf node reached -> evaluate affine function
    disp = $(fname)_hp_list[id]*($(uppercase(fname))_N_PARAMETER+1)*$(uppercase(fname))_N_SOLUTION;
    for(i=0; i < $(uppercase(fname))_N_SOLUTION; i++){
        for(j=0, val=0; j < $(uppercase(fname))_N_PARAMETER; j++)
            val += parameter[j] * $(fname)_feedbacks[disp++];
        val += $(fname)_feedbacks[disp++];
        solution[i] = val;
    }
}
          """)
    close(fsrc)
end
