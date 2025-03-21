:(function (ˍ₋out, x)
      begin
          ˍ₋out_1 = (view)(ˍ₋out, 2:5, 2:5)
          var"##323" = for var"##325" = (zip)(1:4, (Symbolics.reset_to_one)(1:4))
                  begin
                      j = var"##325"[1]
                      j′ = var"##325"[2]
                      for var"##324" = (zip)(1:4, (Symbolics.reset_to_one)(1:4))
                          begin
                              i = var"##324"[1]
                              i′ = var"##324"[2]
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