:(function (ˍ₋out, x)
      begin
          ˍ₋out_1 = (view)(ˍ₋out, 1:5, 1:5)
          var"##294" = (Symbolics.broadcast_assign!)(ˍ₋out_1, 0)
          ˍ₋out_2 = (view)(ˍ₋out, 2:4, 2:4)
          var"##295" = for var"##297" = (zip)(2:4, (Symbolics.reset_to_one)(2:4))
                  begin
                      j = var"##297"[1]
                      j′ = var"##297"[2]
                      for var"##296" = (zip)(2:4, (Symbolics.reset_to_one)(2:4))
                          begin
                              i = var"##296"[1]
                              i′ = var"##296"[2]
                              begin
                                  ˍ₋out_2[(CartesianIndex)(i′, j′)] = (+)((getindex)(ˍ₋out_2, i′, j′), (*)(1//2, (+)((+)((+)((getindex)(x, (+)(-1, i), j), (getindex)(x, (+)(1, i), j)), (getindex)(x, i, (+)(-1, j))), (getindex)(x, i, (+)(1, j)))))
                                  nothing
                              end
                          end
                      end
                  end
              end
          nothing
      end
  end)