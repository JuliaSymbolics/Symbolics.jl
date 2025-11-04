:(function (u,)
      begin
          _out = (zeros)(Float64, (5, 5))
          var"%_out" = for _2 = 1:1:5
                  for _1 = 1:1:5
                      begin
                          _out[(CartesianIndex)(_1, _2)] = (+)((getindex)(_out, _1, _2), (getindex)(u, (Main.limit2)((+)(-1, (getindex)(1:1:5, _2)), 5), (Main.limit2)((+)(1, (getindex)(1:1:5, _1)), 5)))
                          nothing
                      end
                  end
              end
          _out
      end
  end)