
pro get_good_blocks, coords, r_max, tree, GOODBLOCKS=good_blocks, CBLOCKS=cblocks


centerx = coords[0]
centery = coords[1]
centerz = coords[2]

range_xl = centerx - r_max
range_xu = centerx + r_max
range_yl = centery - r_max
range_yu = centery + r_max
range_zl = centerz - r_max
range_zu = centerz + r_max

good_blocks = where(tree.nodetype EQ 1 AND $
                    tree.bndbox[0,0] LE range_xu AND tree.bndbox[1,0] GT range_xl AND $
                    tree.bndbox[0,1] LE range_yu AND tree.bndbox[1,1] GT range_yl AND $
                    tree.bndbox[0,2] LE range_zu AND tree.bndbox[1,2] GT range_zl, cblocks)


end
