## Printing
function print_ws(ws,j)
    print("\r>> #$j\
          |Pending : $(length(ws.Sdown)+length(ws.Sup))\
          |Finished: $(length(ws.F))|       ");
end
function print_final(ws)
    printstyled("\r======= \
           Solution with $(length(ws.F)) critical regions \
          =======     \n",color=:light_magenta);
end

function print_tree_build(U,jump_list,depth)
    print("\r>>\
          |# nodes: $(length(jump_list))\
          |depth: $depth\
          |pending : $(length(U))|         ")
end

function Base.show(io::IO, sol::Solution)
    println("==== Parametric Solution ========")
    if(length(sol.CRs) == 0)
        println(io, "Empty parametric solution")
        println(io, "Status: $(sol.status)")
    else
        println(io, "Status: $(sol.status)")
        println(io, "# of regions:     $(length(sol.CRs))")
        println(io, "# of parameters:  $(sol.problem.n_theta)")
        println(io, "# of constraints: $(length(sol.problem.norm_factors))")
        println(io, "# of variables:   $(sol.problem.n)")
        println(io, "# of output dims: $(size(sol.CRs[1].z,2))")
    end
    println("=================================")
end

function Base.show(io::IO, tree::BinarySearchTree)
    println("====== Binary Search Tree ======")
    println(io, "Tree depth:       $(tree.depth)")
    println(io, "# of Nodes:       $(length(tree.jump_list))")

    n_real = length(tree.halfplanes)
    !isempty(tree.feedbacks) && (n_real += length(tree.feedbacks)*length(tree.feedbacks[1]))
    !isempty(tree.duals) && (n_real+=length(tree.duals)*length(tree.duals[1]))

    println(io, "# of Reals:       $n_real")
    println(io, "# of Integers:    $(2*length(tree.jump_list))")
    println(io, "# of Feedbacks:   $(length(tree.feedbacks))")
    println(io, "# of Half-planes: $(size(tree.halfplanes,2))")

    println("================================")
end
