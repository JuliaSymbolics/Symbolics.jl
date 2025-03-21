:(function (ˍ₋out, u)
      begin
          ˍ₋out_input_1 = let _out = (zeros)(Float64, (map)(length, (Base.OneTo(5), Base.OneTo(5)))), var"##316" = for var"##318" = (zip)(1:5, (Symbolics.reset_to_one)(1:5))
                          begin
                              j = var"##318"[1]
                              j′ = var"##318"[2]
                              for var"##317" = (zip)(1:5, (Symbolics.reset_to_one)(1:5))
                                  begin
                                      i = var"##317"[1]
                                      i′ = var"##317"[2]
                                      begin
                                          _out[(CartesianIndex)(i′, j′)] = (+)((getindex)(_out, i′, j′), (getindex)(u, (Main.limit2)((+)(-1, i), 5), (Main.limit2)((+)(1, j), 5)))
                                          nothing
                                      end
                                  end
                              end
                          end
                      end
                  _out
              end
          for var"##309" = (zip)(1:5, (Symbolics.reset_to_one)(1:5))
              begin
                  j = var"##309"[1]
                  j′ = var"##309"[2]
                  for var"##308" = (zip)(1:5, (Symbolics.reset_to_one)(1:5))
                      begin
                          i = var"##308"[1]
                          i′ = var"##308"[2]
                          begin
                              ˍ₋out[(CartesianIndex)(i′, j′)] = (+)((getindex)(ˍ₋out, i′, j′), (getindex)(ˍ₋out_input_1, j, i))
                              nothing
                          end
                      end
                  end
              end
          end
      end
  end)