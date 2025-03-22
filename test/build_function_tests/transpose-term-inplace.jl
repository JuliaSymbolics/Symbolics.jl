:(function (ˍ₋out, x)
      for var"%jj′" = (zip)(1:4, (Symbolics.reset_to_one)(1:4))
          begin
              j = var"%jj′"[1]
              j′ = var"%jj′"[2]
              for var"%ii′" = (zip)(1:4, (Symbolics.reset_to_one)(1:4))
                  begin
                      i = var"%ii′"[1]
                      i′ = var"%ii′"[2]
                      begin
                          ˍ₋out[(CartesianIndex)(i′, j′)] = (+)((getindex)(ˍ₋out, i′, j′), (getindex)(x, j, i))
                          nothing
                      end
                  end
              end
          end
      end
  end)