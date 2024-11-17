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
