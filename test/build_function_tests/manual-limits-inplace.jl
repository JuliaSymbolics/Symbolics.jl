:(function (ˍ₋out, u)
      for var"##349" = (zip)(1:5, (Symbolics.reset_to_one)(1:5))
          begin
              j = var"##349"[1]
              j′ = var"##349"[2]
              for var"##348" = (zip)(1:5, (Symbolics.reset_to_one)(1:5))
                  begin
                      i = var"##348"[1]
                      i′ = var"##348"[2]
                      begin
                          ˍ₋out[(CartesianIndex)(i′, j′)] = (+)((getindex)(ˍ₋out, i′, j′), (getindex)(u, (Main.limit2)((+)(-1, i), 5), (Main.limit2)((+)(1, j), 5)))
                          nothing
                      end
                  end
              end
          end
      end
  end)