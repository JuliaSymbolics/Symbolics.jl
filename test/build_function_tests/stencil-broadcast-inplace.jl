:(function (ˍ₋out, x)
      begin
          ˍ₋out_2_input_1 = (broadcast)(+, x, (adjoint)(x))
          ˍ₋out_1 = (view)(ˍ₋out, 1:6, 1:6)
          var"%ˍ₋out_1" = (Symbolics.broadcast_assign!)(ˍ₋out_1, 0)
          ˍ₋out_2 = (view)(ˍ₋out, 2:5, 2:5)
          var"%ˍ₋out_2" = for var"%jj′" = (zip)(1:4, (Symbolics.reset_to_one)(1:4))
                  begin
                      j = var"%jj′"[1]
                      j′ = var"%jj′"[2]
                      for var"%ii′" = (zip)(1:4, (Symbolics.reset_to_one)(1:4))
                          begin
                              i = var"%ii′"[1]
                              i′ = var"%ii′"[2]
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