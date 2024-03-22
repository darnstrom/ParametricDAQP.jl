## Printing
function print_ws(ws,j)
    print("\r>> #$j\
          |Down: $(length(ws.Sdown))\
          |Up: $(length(ws.Sup))\
          |Fin: $(length(ws.F))|       ");
end
function print_final(ws)
    print("\r======= |\
          | Fin: $(length(ws.F)) \
          | # LPs: $(ws.nLPs) \
          | explored: $(length(ws.explored)) \
          ||=======     \n");
end
