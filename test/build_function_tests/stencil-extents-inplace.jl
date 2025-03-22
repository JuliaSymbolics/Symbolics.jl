:(function (ˍ₋out, x)
      begin
          ˍ₋out_1 = (view)(ˍ₋out, 1:5, 1:5)
          var"%ˍ₋out_1" = (Symbolics.broadcast_assign!)(ˍ₋out_1, 0)
          ˍ₋out_2 = (view)(ˍ₋out, 2:4, 2:4)
          var"%ˍ₋out_2" = for var"%jj′" = (zip)(2:4, (Symbolics.reset_to_one)(2:4))
                  begin
                      j = var"%jj′"[1]
                      j′ = var"%jj′"[2]
                      for var"%ii′" = (zip)(2:4, (Symbolics.reset_to_one)(2:4))
                          begin
                              i = var"%ii′"[1]
                              i′ = var"%ii′"[2]
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