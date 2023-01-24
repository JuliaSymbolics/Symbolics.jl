:(function (x,)
      let ˍ₋out = zeros(Float64, map(length, (1:5, 1:5)))
          begin
              ˍ₋out_1 = (view)(ˍ₋out, 1:5, 1:5)
              ˍ₋out_1 .= 0
              ˍ₋out_2 = (view)(ˍ₋out, 2:4, 2:4)
              for (j, j′) in zip(2:4, reset_to_one(2:4))
                  for (i, i′) in zip(2:4, reset_to_one(2:4))
                      ˍ₋out_2[i′, j′] = (+)(ˍ₋out_2[i′, j′],
                                            (*)(1 // 2,
                                                (+)((+)((+)((getindex)(x, i, (+)(1, j)),
                                                            (getindex)(x, i, (+)(-1, j))),
                                                        (getindex)(x, (+)(-1, i), j)),
                                                    (getindex)(x, (+)(1, i), j))))
                  end
              end
          end
          ˍ₋out
      end
  end)
