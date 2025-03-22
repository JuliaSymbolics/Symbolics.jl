:(function (x,)
      let _out = (zeros)(Float64, (map)(length, (1:5, 1:5))), var"%_out" = begin
                  _out_1 = (view)(_out, 1:5, 1:5)
                  var"%_out_1" = (Symbolics.broadcast_assign!)(_out_1, 0)
                  _out_2 = (view)(_out, 2:4, 2:4)
                  var"%_out_2" = for var"%jj′" = (zip)(2:4, (Symbolics.reset_to_one)(2:4))
                          begin
                              j = var"%jj′"[1]
                              j′ = var"%jj′"[2]
                              for var"%ii′" = (zip)(2:4, (Symbolics.reset_to_one)(2:4))
                                  begin
                                      i = var"%ii′"[1]
                                      i′ = var"%ii′"[2]
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