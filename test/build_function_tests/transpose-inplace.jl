:(function (ˍ₋out, x)
      for var"##317" = (zip)(1:4, (Symbolics.reset_to_one)(1:4))
          begin
              j = var"##317"[1]
              j′ = var"##317"[2]
              for var"##316" = (zip)(1:4, (Symbolics.reset_to_one)(1:4))
                  begin
                      i = var"##316"[1]
                      i′ = var"##316"[2]
                      begin
                          ˍ₋out[(CartesianIndex)(i′, j′)] = (+)((getindex)(ˍ₋out, i′, j′), (getindex)(x, j, i))
                          nothing
                      end
                  end
              end
          end
      end
  end)