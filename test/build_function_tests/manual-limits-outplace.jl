:(function (ˍ₋out, u)
      begin
          for (j, j′) = zip(1:5, reset_to_one(1:5))
              for (i, i′) = zip(1:5, reset_to_one(1:5))
                  ˍ₋out[i′, j′] = (+)(ˍ₋out[i′, j′], (getindex)(u, (limit)((+)(-1, i), 5), (limit)((+)(1, j), 5)))
              end
          end
      end
  end)