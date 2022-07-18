:(function (ˍ₋out, u)
      begin
          ˍ₋out_input_1 = let _out = zeros(Float64, map(length, (Base.OneTo(5), Base.OneTo(5))))
                  begin
                      for (j, j′) = zip(1:5, reset_to_one(1:5))
                          for (i, i′) = zip(1:5, reset_to_one(1:5))
                              _out[i′, j′] = (+)(_out[i′, j′], (getindex)(u, (limit)((+)(-1, i), 5), (limit)((+)(1, j), 5)))
                          end
                      end
                  end
                  _out
              end
          for (j, j′) = zip(Base.OneTo(5), reset_to_one(Base.OneTo(5)))
              for (i, i′) = zip(Base.OneTo(5), reset_to_one(Base.OneTo(5)))
                  ˍ₋out[i′, j′] = (+)(ˍ₋out[i′, j′], (getindex)(ˍ₋out_input_1, j, i))
              end
          end
      end
  end)