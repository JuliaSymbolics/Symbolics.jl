:(function (ˍ₋out, x)
      begin
          for (j, j′) = zip(1:4, reset_to_one(1:4))
              for (i, i′) = zip(1:4, reset_to_one(1:4))
                  ˍ₋out[i′, j′] = (+)(ˍ₋out[i′, j′], (getindex)(x, j, i))
              end
          end
      end
  end)