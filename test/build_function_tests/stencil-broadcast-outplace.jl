:(function (x,)
      let _out = (zeros)(Float64, (map)(length, (1:6, 1:6))), var"##289" = begin
                  _out_2_input_1 = (broadcast)(+, x, (adjoint)(x))
                  _out_1 = (view)(_out, 1:6, 1:6)
                  var"##290" = (Symbolics.broadcast_assign!)(_out_1, 0)
                  _out_2 = (view)(_out, 2:5, 2:5)
                  var"##291" = for var"##293" = (zip)(Base.OneTo(4), (Symbolics.reset_to_one)(Base.OneTo(4)))
                          begin
                              j = var"##293"[1]
                              j′ = var"##293"[2]
                              for var"##292" = (zip)(Base.OneTo(4), (Symbolics.reset_to_one)(Base.OneTo(4)))
                                  begin
                                      i = var"##292"[1]
                                      i′ = var"##292"[2]
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