"""A single objective for optimization with Krotov's method."""
struct Objective
    intial_state
    H
    target
end
