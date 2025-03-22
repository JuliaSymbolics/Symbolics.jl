:(function (x,)
      let _out = (zeros)(Float64, (map)(length, (1:6, 1:6))), var"%_out" = begin
                  _out_2_input_1 = (broadcast)(+, x, (adjoint)(x))
                  _out_1 = (view)(_out, 1:6, 1:6)
                  var"%_out_1" = (Symbolics.broadcast_assign!)(_out_1, 0)
                  _out_2 = (view)(_out, 2:5, 2:5)
                  var"%_out_2" = for var"%jj′" = (zip)(Base.OneTo(4), (Symbolics.reset_to_one)(Base.OneTo(4)))
                          begin
                              j = var"%jj′"[1]
                              j′ = var"%jj′"[2]
                              for var"%ii′" = (zip)(Base.OneTo(4), (Symbolics.reset_to_one)(Base.OneTo(4)))
                                  begin
                                      i = var"%ii′"[1]
                                      i′ = var"%ii′"[2]
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