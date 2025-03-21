:(function (x,)
      let _out = (zeros)(Float64, (map)(length, (Base.OneTo(4), Base.OneTo(4)))), var"##273" = for var"##275" = (zip)(1:4, (Symbolics.reset_to_one)(1:4))
                  begin
                      j = var"##275"[1]
                      j′ = var"##275"[2]
                      for var"##274" = (zip)(1:4, (Symbolics.reset_to_one)(1:4))
                          begin
                              i = var"##274"[1]
                              i′ = var"##274"[2]
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