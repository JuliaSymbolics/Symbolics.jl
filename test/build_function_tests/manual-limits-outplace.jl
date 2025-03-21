:(function (u,)
      let _out = (zeros)(Float64, (map)(length, (Base.OneTo(5), Base.OneTo(5)))), var"##350" = for var"##352" = (zip)(1:5, (Symbolics.reset_to_one)(1:5))
                  begin
                      j = var"##352"[1]
                      j′ = var"##352"[2]
                      for var"##351" = (zip)(1:5, (Symbolics.reset_to_one)(1:5))
                          begin
                              i = var"##351"[1]
                              i′ = var"##351"[2]
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
  end)