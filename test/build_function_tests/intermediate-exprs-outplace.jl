:(function (u,)
      let _out = (zeros)(Float64, (map)(length, (Base.OneTo(5), Base.OneTo(5)))), var"##355" = begin
                  _out_input_1 = let _out = (zeros)(Float64, (map)(length, (Base.OneTo(5), Base.OneTo(5)))), var"##358" = for var"##360" = (zip)(1:5, (Symbolics.reset_to_one)(1:5))
                                  begin
                                      j = var"##360"[1]
                                      j′ = var"##360"[2]
                                      for var"##359" = (zip)(1:5, (Symbolics.reset_to_one)(1:5))
                                          begin
                                              i = var"##359"[1]
                                              i′ = var"##359"[2]
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
                  for var"##357" = (zip)(1:5, (Symbolics.reset_to_one)(1:5))
                      begin
                          j = var"##357"[1]
                          j′ = var"##357"[2]
                          for var"##356" = (zip)(1:5, (Symbolics.reset_to_one)(1:5))
                              begin
                                  i = var"##356"[1]
                                  i′ = var"##356"[2]
                                  begin
                                      _out[(CartesianIndex)(i′, j′)] = (+)((getindex)(_out, i′, j′), (getindex)(_out_input_1, j, i))
                                      nothing
                                  end
                              end
                          end
                      end
                  end
              end
          _out
      end
  end)