:(function (x,)
      let _out = (zeros)(Float64, (map)(length, (Base.OneTo(4), Base.OneTo(4)))), var"##318" = for var"##320" = (zip)(1:4, (Symbolics.reset_to_one)(1:4))
                  begin
                      j = var"##320"[1]
                      j′ = var"##320"[2]
                      for var"##319" = (zip)(1:4, (Symbolics.reset_to_one)(1:4))
                          begin
                              i = var"##319"[1]
                              i′ = var"##319"[2]
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