:(function (x,)
      let _out = (zeros)(Float64, (map)(length, (1:6, 1:6))), var"##334" = begin
                  _out_2_input_1 = (broadcast)(+, x, (adjoint)(x))
                  _out_1 = (view)(_out, 1:6, 1:6)
                  var"##335" = (Symbolics.broadcast_assign!)(_out_1, 0)
                  _out_2 = (view)(_out, 2:5, 2:5)
                  var"##336" = for var"##338" = (zip)(Base.OneTo(4), (Symbolics.reset_to_one)(Base.OneTo(4)))
                          begin
                              j = var"##338"[1]
                              j′ = var"##338"[2]
                              for var"##337" = (zip)(Base.OneTo(4), (Symbolics.reset_to_one)(Base.OneTo(4)))
                                  begin
                                      i = var"##337"[1]
                                      i′ = var"##337"[2]
                                      begin
                                          _out_2[(CartesianIndex)(i′, j′)] = (+)((getindex)(_out_2, i′, j′), (+)(1, (getindex)(_out_2_input_1, i, j)))
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