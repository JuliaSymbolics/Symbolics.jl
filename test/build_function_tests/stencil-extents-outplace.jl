:(function (x,)
      let _out = (zeros)(Float64, (map)(length, (1:5, 1:5))), var"##298" = begin
                  _out_1 = (view)(_out, 1:5, 1:5)
                  var"##299" = (Symbolics.broadcast_assign!)(_out_1, 0)
                  _out_2 = (view)(_out, 2:4, 2:4)
                  var"##300" = for var"##302" = (zip)(2:4, (Symbolics.reset_to_one)(2:4))
                          begin
                              j = var"##302"[1]
                              j′ = var"##302"[2]
                              for var"##301" = (zip)(2:4, (Symbolics.reset_to_one)(2:4))
                                  begin
                                      i = var"##301"[1]
                                      i′ = var"##301"[2]
                                      begin
                                          _out_2[(CartesianIndex)(i′, j′)] = (+)((getindex)(_out_2, i′, j′), (*)(1//2, (+)((+)((+)((getindex)(x, (+)(-1, i), j), (getindex)(x, (+)(1, i), j)), (getindex)(x, i, (+)(-1, j))), (getindex)(x, i, (+)(1, j)))))
                                          nothing
                                      end
                                  end
                              end
                          end
                      end
                  nothing
              end
          _out
      end
  end)