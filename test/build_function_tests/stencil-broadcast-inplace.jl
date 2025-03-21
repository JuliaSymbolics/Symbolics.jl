:(function (ˍ₋out, x)
      begin
          ˍ₋out_2_input_1 = (broadcast)(+, x, (adjoint)(x))
          ˍ₋out_1 = (view)(ˍ₋out, 1:6, 1:6)
          var"##285" = (Symbolics.broadcast_assign!)(ˍ₋out_1, 0)
          ˍ₋out_2 = (view)(ˍ₋out, 2:5, 2:5)
          var"##286" = for var"##288" = (zip)(Base.OneTo(4), (Symbolics.reset_to_one)(Base.OneTo(4)))
                  begin
                      j = var"##288"[1]
                      j′ = var"##288"[2]
                      for var"##287" = (zip)(Base.OneTo(4), (Symbolics.reset_to_one)(Base.OneTo(4)))
                          begin
                              i = var"##287"[1]
                              i′ = var"##287"[2]
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