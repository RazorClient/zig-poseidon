const std = @import("std");

test "babyBear16" {
    std.testing.log_level = .debug;
     std.debug.print("***Main.zig test***",.{});
    _ = @import("instances/babybear16.zig");
}
