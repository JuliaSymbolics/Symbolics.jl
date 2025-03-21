:(function (ˍ₋out, x)
      begin
          ˍ₋out_1 = (view)(ˍ₋out, 2:5, 2:5)
          var"##278" = for var"##280" = (zip)(1:4, (Symbolics.reset_to_one)(1:4))
                  begin
                      j = var"##280"[1]
                      j′ = var"##280"[2]
                      for var"##279" = (zip)(1:4, (Symbolics.reset_to_one)(1:4))
                          begin
                              i = var"##279"[1]
                              i′ = var"##279"[2]
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