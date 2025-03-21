:(function (ˍ₋out, u)
      begin
          ˍ₋out_input_1 = let _out = (zeros)(Float64, (map)(length, (Base.OneTo(5), Base.OneTo(5)))), var"##361" = for var"##363" = (zip)(1:5, (Symbolics.reset_to_one)(1:5))
                          begin
                              j = var"##363"[1]
                              j′ = var"##363"[2]
                              for var"##362" = (zip)(1:5, (Symbolics.reset_to_one)(1:5))
                                  begin
                                      i = var"##362"[1]
                                      i′ = var"##362"[2]
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
          for var"##354" = (zip)(1:5, (Symbolics.reset_to_one)(1:5))
              begin
                  j = var"##354"[1]
                  j′ = var"##354"[2]
                  for var"##353" = (zip)(1:5, (Symbolics.reset_to_one)(1:5))
                      begin
                          i = var"##353"[1]
                          i′ = var"##353"[2]
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