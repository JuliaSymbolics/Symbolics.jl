:(function (ˍ₋out, x)
      begin
          ˍ₋out_2_input_1 = (broadcast)(+, x, (adjoint)(x))
          ˍ₋out_1 = (view)(ˍ₋out, 1:6, 1:6)
          var"##330" = (Symbolics.broadcast_assign!)(ˍ₋out_1, 0)
          ˍ₋out_2 = (view)(ˍ₋out, 2:5, 2:5)
          var"##331" = for var"##333" = (zip)(Base.OneTo(4), (Symbolics.reset_to_one)(Base.OneTo(4)))
                  begin
                      j = var"##333"[1]
                      j′ = var"##333"[2]
                      for var"##332" = (zip)(Base.OneTo(4), (Symbolics.reset_to_one)(Base.OneTo(4)))
                          begin
                              i = var"##332"[1]
                              i′ = var"##332"[2]
                              begin
                                  ˍ₋out_2[(CartesianIndex)(i′, j′)] = (+)((getindex)(ˍ₋out_2, i′, j′), (+)(1, (getindex)(ˍ₋out_2_input_1, i, j)))
                                  nothing
                              end
                          end
                      end
                  end
              end
          nothing
      end
  end)