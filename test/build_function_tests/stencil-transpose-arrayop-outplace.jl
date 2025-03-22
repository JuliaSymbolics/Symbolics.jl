:(function (x,)
      let _out = (zeros)(Float64, (map)(length, (1:6, 1:6))), var"%_out" = begin
                  _out_1 = (view)(_out, 2:5, 2:5)
                  var"%_out_1" = for var"%jj′" = (zip)(1:4, (Symbolics.reset_to_one)(1:4))
                          begin
                              j = var"%jj′"[1]
                              j′ = var"%jj′"[2]
                              for var"%ii′" = (zip)(1:4, (Symbolics.reset_to_one)(1:4))
                                  begin
                                      i = var"%ii′"[1]
                                      i′ = var"%ii′"[2]
                                      begin
                                          _out_1[(CartesianIndex)(i′, j′)] = (+)((getindex)(_out_1, i′, j′), (getindex)(x, j, i))
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