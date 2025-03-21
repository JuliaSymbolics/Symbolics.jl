:(function (u,)
      let _out = (zeros)(Float64, (map)(length, (Base.OneTo(5), Base.OneTo(5)))), var"##310" = begin
                  _out_input_1 = let _out = (zeros)(Float64, (map)(length, (Base.OneTo(5), Base.OneTo(5)))), var"##313" = for var"##315" = (zip)(1:5, (Symbolics.reset_to_one)(1:5))
                                  begin
                                      j = var"##315"[1]
                                      j′ = var"##315"[2]
                                      for var"##314" = (zip)(1:5, (Symbolics.reset_to_one)(1:5))
                                          begin
                                              i = var"##314"[1]
                                              i′ = var"##314"[2]
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
                  for var"##312" = (zip)(1:5, (Symbolics.reset_to_one)(1:5))
                      begin
                          j = var"##312"[1]
                          j′ = var"##312"[2]
                          for var"##311" = (zip)(1:5, (Symbolics.reset_to_one)(1:5))
                              begin
                                  i = var"##311"[1]
                                  i′ = var"##311"[2]
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