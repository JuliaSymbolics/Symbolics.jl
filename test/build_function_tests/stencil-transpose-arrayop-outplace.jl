:(function (x,)
      let _out = (zeros)(Float64, (map)(length, (1:6, 1:6))), var"##281" = begin
                  _out_1 = (view)(_out, 2:5, 2:5)
                  var"##282" = for var"##284" = (zip)(1:4, (Symbolics.reset_to_one)(1:4))
                          begin
                              j = var"##284"[1]
                              j′ = var"##284"[2]
                              for var"##283" = (zip)(1:4, (Symbolics.reset_to_one)(1:4))
                                  begin
                                      i = var"##283"[1]
                                      i′ = var"##283"[2]
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