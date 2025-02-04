:(function (ˍ₋out, u)
      begin
          for (j, j′) = zip(1:5, reset_to_one(1:5))
              for (i, i′) = zip(1:5, reset_to_one(1:5))
                  ˍ₋out[i′, j′] = (+)(ˍ₋out[i′, j′], (getindex)(u, (Main.limit2)((+)(-1, i), 5), (Main.limit2)((+)(1, j), 5)))
              end
          end
      end
  end)