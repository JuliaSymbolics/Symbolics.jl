:(function (ˍ₋out, u)
      begin
          ˍ₋out_input_1 = begin
                  _out = (zeros)(Float64, (map)(length, (Base.OneTo(5), Base.OneTo(5))))
                  var"%_out" = for var"%jj′" = (zip)(Base.OneTo(5), (Symbolics.reset_to_one)(Base.OneTo(5)))
                          begin
                              j = var"%jj′"[1]
                              j′ = var"%jj′"[2]
                              for var"%ii′" = (zip)(Base.OneTo(5), (Symbolics.reset_to_one)(Base.OneTo(5)))
                                  begin
                                      i = var"%ii′"[1]
                                      i′ = var"%ii′"[2]
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
          for var"%jj′" = (zip)(Base.OneTo(5), (Symbolics.reset_to_one)(Base.OneTo(5)))
              begin
                  j = var"%jj′"[1]
                  j′ = var"%jj′"[2]
                  for var"%ii′" = (zip)(Base.OneTo(5), (Symbolics.reset_to_one)(Base.OneTo(5)))
                      begin
                          i = var"%ii′"[1]
                          i′ = var"%ii′"[2]
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