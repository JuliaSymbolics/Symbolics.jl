:(function (ˍ₋out, x)
      begin for (j, j′) in zip(1:4, reset_to_one(1:4))
          for (i, i′) in zip(1:4, reset_to_one(1:4))
              ˍ₋out[i′, j′] = (+)(ˍ₋out[i′, j′], (getindex)(x, j, i))
          end
      end end
  end)
