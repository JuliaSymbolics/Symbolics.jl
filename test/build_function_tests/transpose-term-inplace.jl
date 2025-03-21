:(function (ˍ₋out, x)
      for var"##277" = (zip)(1:4, (Symbolics.reset_to_one)(1:4))
          begin
              j = var"##277"[1]
              j′ = var"##277"[2]
              for var"##276" = (zip)(1:4, (Symbolics.reset_to_one)(1:4))
                  begin
                      i = var"##276"[1]
                      i′ = var"##276"[2]
                      begin
                          ˍ₋out[(CartesianIndex)(i′, j′)] = (+)((getindex)(ˍ₋out, i′, j′), (getindex)(x, j, i))
                          nothing
                      end
                  end
              end
          end
      end
  end)