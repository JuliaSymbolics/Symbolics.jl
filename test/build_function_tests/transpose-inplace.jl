:(function (ˍ₋out, x)
      begin
          _out = ˍ₋out
          var"%_out" = for _2 = 1:1:4
                  for _1 = 1:1:4
                      begin
                          _out[(CartesianIndex)(_1, _2)] = (+)((getindex)(_out, _1, _2), (getindex)(x, _2, _1))
                          nothing
                      end
                  end
              end
          _out
      end
  end)