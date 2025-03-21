:(function (x,)
      let _out = (zeros)(Float64, (map)(length, (1:6, 1:6))), var"##326" = begin
                  _out_1 = (view)(_out, 2:5, 2:5)
                  var"##327" = for var"##329" = (zip)(1:4, (Symbolics.reset_to_one)(1:4))
                          begin
                              j = var"##329"[1]
                              j′ = var"##329"[2]
                              for var"##328" = (zip)(1:4, (Symbolics.reset_to_one)(1:4))
                                  begin
                                      i = var"##328"[1]
                                      i′ = var"##328"[2]
                                      begin
                                          _out_1[(CartesianIndex)(i′, j′)] = (+)((getindex)(_out_1, i′, j′), (getindex)(x, j, i))
                                          nothing
                                      end
                                  end
                              end
                          end
                      end
                  nothing
              end
          _out
      end
  end)