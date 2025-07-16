:(function (ˍ₋out, x)
      begin
          ˍ₋out_1 = (view)(ˍ₋out, 2:5, 2:5)
          var"%ˍ₋out_1" = for var"%jj′" = (zip)(Base.OneTo(4), (Symbolics.reset_to_one)(Base.OneTo(4)))
                  begin
                      j = var"%jj′"[1]
                      j′ = var"%jj′"[2]
                      for var"%ii′" = (zip)(Base.OneTo(4), (Symbolics.reset_to_one)(Base.OneTo(4)))
                          begin
                              i = var"%ii′"[1]
                              i′ = var"%ii′"[2]
                              begin
                                  ˍ₋out_1[(CartesianIndex)(i′, j′)] = (+)((getindex)(ˍ₋out_1, i′, j′), (getindex)(x, j, i))
                                  nothing
                              end
                          end
                      end
                  end
              end
          nothing
      end
  end)