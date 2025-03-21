:(function (ˍ₋out, x)
      for var"##272" = (zip)(1:4, (Symbolics.reset_to_one)(1:4))
          begin
              j = var"##272"[1]
              j′ = var"##272"[2]
              for var"##271" = (zip)(1:4, (Symbolics.reset_to_one)(1:4))
                  begin
                      i = var"##271"[1]
                      i′ = var"##271"[2]
                      begin
                          ˍ₋out[(CartesianIndex)(i′, j′)] = (+)((getindex)(ˍ₋out, i′, j′), (getindex)(x, j, i))
                          nothing
                      end
                  end
              end
          end
      end
  end)