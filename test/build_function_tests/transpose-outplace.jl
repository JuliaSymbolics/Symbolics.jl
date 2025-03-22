:(function (x,)
      let _out = (zeros)(Float64, (map)(length, (Base.OneTo(4), Base.OneTo(4)))), var"%_out" = for var"%jj′" = (zip)(1:4, (Symbolics.reset_to_one)(1:4))
                  begin
                      j = var"%jj′"[1]
                      j′ = var"%jj′"[2]
                      for var"%ii′" = (zip)(1:4, (Symbolics.reset_to_one)(1:4))
                          begin
                              i = var"%ii′"[1]
                              i′ = var"%ii′"[2]
                              begin
                                  _out[(CartesianIndex)(i′, j′)] = (+)((getindex)(_out, i′, j′), (getindex)(x, j, i))
                                  nothing
                              end
                          end
                      end
                  end
              end
          _out
      end
  end)