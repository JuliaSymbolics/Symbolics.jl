:(function (ˍ₋out, x)
      begin
          ˍ₋out_1 = (view)(ˍ₋out, 1:5, 1:5)
          ˍ₋out_1 .= 0
          ˍ₋out_2 = (view)(ˍ₋out, 2:4, 2:4)
          for (j, j′) = zip(2:4, reset_to_one(2:4))
              for (i, i′) = zip(2:4, reset_to_one(2:4))
                  ˍ₋out_2[i′, j′] = (+)(ˍ₋out_2[i′, j′], (/)((+)((+)((+)((getindex)(x, i, (+)(-1, j)), (getindex)(x, i, (+)(1, j))), (getindex)(x, (+)(-1, i), j)), (getindex)(x, (+)(1, i), j)), 2))
              end
          end
      end
  end)