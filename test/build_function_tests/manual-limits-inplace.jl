:(function (ˍ₋out, u)
      for var"%jj′" = (zip)(Base.OneTo(5), (Symbolics.reset_to_one)(Base.OneTo(5)))
          begin
              j = var"%jj′"[1]
              j′ = var"%jj′"[2]
              for var"%ii′" = (zip)(Base.OneTo(5), (Symbolics.reset_to_one)(Base.OneTo(5)))
                  begin
                      i = var"%ii′"[1]
                      i′ = var"%ii′"[2]
                      begin
                          ˍ₋out[(CartesianIndex)(i′, j′)] = (+)((getindex)(ˍ₋out, i′, j′), (getindex)(u, (Main.limit2)((+)(-1, i), 5), (Main.limit2)((+)(1, j), 5)))
                          nothing
                      end
                  end
              end
          end
      end
  end)