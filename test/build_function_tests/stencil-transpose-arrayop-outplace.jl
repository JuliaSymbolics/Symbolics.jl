:(function (ˍ₋out, x)
      begin
          ˍ₋out_1 = (view)(ˍ₋out, 2:5, 2:5)
          for (j, j′) = zip(1:4, reset_to_one(1:4))
              for (i, i′) = zip(1:4, reset_to_one(1:4))
                  ˍ₋out_1[i′, j′] = (+)(ˍ₋out_1[i′, j′], (getindex)(x, j, i))
              end
          end
      end
  end)