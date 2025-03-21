:(function (ˍ₋out, x)
      for var"##322" = (zip)(1:4, (Symbolics.reset_to_one)(1:4))
          begin
              j = var"##322"[1]
              j′ = var"##322"[2]
              for var"##321" = (zip)(1:4, (Symbolics.reset_to_one)(1:4))
                  begin
                      i = var"##321"[1]
                      i′ = var"##321"[2]
                      begin
                          ˍ₋out[(CartesianIndex)(i′, j′)] = (+)((getindex)(ˍ₋out, i′, j′), (getindex)(x, j, i))
                          nothing
                      end
                  end
              end
          end
      end
  end)