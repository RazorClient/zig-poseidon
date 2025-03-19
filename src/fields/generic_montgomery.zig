const std = @import("std");

// MontgomeryField31 is a generic implementation of a field with modulus with at most 31 bits.
pub fn MontgomeryField31(comptime modulus: u32) type {
    const R: u64 = 1 << 32;
    const R_square_mod_modulus: u64 = @intCast((@as(u128, R) * @as(u128, R)) % modulus);

    // modulus_prime = -modulus^-1 mod R
    const modulus_prime = R - euclideanAlgorithm(modulus, R) % R;
    std.debug.assert(modulus * modulus_prime % R == R - 1);

    return struct {
        const Self = @This();
        pub const PrimeModulus = modulus; 
        pub const FieldElem = u32;
        pub const MontFieldElem = struct {
            value: u32,
        };
        pub fn getModulus() u32 {
            return PrimeModulus;
        }
        pub fn toMontgomery(out: *MontFieldElem, value: FieldElem) void {
            out.* = .{ .value = montReduce(@as(u64, value) * R_square_mod_modulus) };
        }

        pub fn square(out1: *MontFieldElem, value: MontFieldElem) void {
            mul(out1, value, value);
        }

        pub fn mul(out1: *MontFieldElem, value: MontFieldElem, arg2: MontFieldElem) void {
            out1.* = .{ .value = montReduce(@as(u64, value.value) * @as(u64, arg2.value)) };
        }

        pub fn add(out1: *MontFieldElem, value: MontFieldElem, arg2: MontFieldElem) void {
            var tmp = value.value + arg2.value;
            if (tmp > modulus) {
                tmp -= modulus;
            }
            out1.* = .{ .value = tmp };
        }

        /// Subtract two Montgomery field elements: out = a - b.
        pub fn sub(out: *MontFieldElem, a: MontFieldElem, b: MontFieldElem) void {
            if (a.value >= b.value) {
                out.* = .{ .value = a.value - b.value };
            } else {
                out.* = .{ .value = a.value + modulus - b.value };
            }
        }

        // Negate a Montgomery field element.
        pub fn neg(out: *MontFieldElem, a: MontFieldElem) void {
            if (a.value == 0) {
                out.* = .{ .value = 0 };
            } else {
                out.* = .{ .value = modulus - a.value };
            }
        }

        pub fn toNormal(self: MontFieldElem) FieldElem {
            return montReduce(@as(u64, self.value));
        }

        fn montReduce(mont_value: u64) FieldElem {
            const tmp = mont_value + (((mont_value & 0xFFFFFFFF) * modulus_prime) & 0xFFFFFFFF) * modulus;
            std.debug.assert(tmp % R == 0);
            const t = tmp >> 32;
            if (t >= modulus) {
                return @intCast(t - modulus);
            }
            return @intCast(t);
        }

        /// Fast exponentiation: out = base^exponent (in Montgomery form).
        pub fn pow(out: *MontFieldElem, base: MontFieldElem, exponent: u32) void {
            var result = Self.fromU32(1);
            var tmp_base = base;
            var exp = exponent;

            while (exp != 0) : (exp /= 2) {
                if ((exp & 1) != 0) {
                    var prod: MontFieldElem = undefined;
                    Self.mul(&prod, result, tmp_base);
                    result = prod;
                }
                if (exp > 1) {
                    var sq: MontFieldElem = undefined;
                    Self.square(&sq, tmp_base);
                    tmp_base = sq;
                }
            }
            out.* = result;
        }

        /// Uses Fermat's Little Theorem: a^(modulus - 2).
        pub fn inverse(out: *MontFieldElem, a: MontFieldElem) void {
            Self.pow(out, a, modulus - 2);
        }

        pub fn fromU32(x: u32) MontFieldElem {
            const reduced = x % modulus;
            var fe: MontFieldElem = undefined;
            Self.toMontgomery(&fe, reduced);
            return fe;
        }

        pub fn fromU8(x: u8) MontFieldElem {
            return Self.fromU32(@as(u32, x));
        }

        pub fn fromBytes(bytes: []const u8) MontFieldElem {
            var value: FieldElem = 0;
            
            // 4 bytes (32 bits) since our field is 31-bit
            const len = @min(bytes.len, 4);
            for (0..len) |i| {
                value |= @as(FieldElem, bytes[i]) << @intCast(8 * i);
            }
            
            value %= modulus;
            
            var result: MontFieldElem = undefined;
            toMontgomery(&result, value);
            return result;
            }
        
        /// Convert a Montgomery field element to a 4-byte array (little-endian).
        pub fn toBytes(self: MontFieldElem, out: *[4]u8) void {
            const normal = Self.toNormal(self);
            out[0] = @intCast(normal & 0xFF);
            out[1] = @intCast((normal >> 8) & 0xFF);
            out[2] = @intCast((normal >> 16) & 0xFF);
            out[3] = @intCast((normal >> 24) & 0xFF);
        }
    };
}

fn euclideanAlgorithm(a: u64, b: u64) u64 {
    var t: i64 = 0;
    var new_t: i64 = 1;
    var r: i64 = @intCast(b);
    var new_r: i64 = @intCast(a);

    while (new_r != 0) {
        const quotient = r / new_r;

        const temp_t = t;
        t = new_t;
        new_t = temp_t - quotient * new_t;

        const temp_r = r;
        r = new_r;
        new_r = temp_r - quotient * new_r;
    }

    if (r != 1) {
        @compileError("modular inverse does not exist");
    }

    if (t < 0) {
        t += @intCast(b);
    }
    return @intCast(t);
}

pub fn equals(M: anytype, a: M.MontFieldElem, b: M.MontFieldElem) bool {
    return a.value == b.value;
}

pub fn zero(M: anytype) M.MontFieldElem {
    return .{ .value = 0 };
}

pub fn one(M: anytype) M.MontFieldElem {
    return M.fromU32(1);
}

pub fn isZero(M: anytype, a: M.MontFieldElem) bool {
    return a.value == 0;
}


test "MontgomeryField31 basic arithmetic tests" {
    // Use a prime modulus: 2^31 - 1 (a Mersenne prime).
    const Modulus = 0x7fffffff; // 2147483647
    const M = MontgomeryField31(Modulus);

    // Test conversion from u32 and utility functions.
    const one_elem = one(M);
    const zero_elem = zero(M);
    std.debug.print("One (normal): {}\n", .{ M.toNormal(one_elem) });
    std.debug.print("Zero (normal): {}\n", .{ M.toNormal(zero_elem) });
    try std.testing.expect(equals(M, one_elem, M.fromU32(1)));
    try std.testing.expect(equals(M, zero_elem, M.fromU32(0)));

    // Test addition: 1 + 1 = 2.
    var sum: M.MontFieldElem = undefined;
    M.add(&sum, one_elem, one_elem);
    std.debug.print("1 + 1 = {}\n", .{ M.toNormal(sum) });
    try std.testing.expect(equals(M, sum, M.fromU32(2)));

    // Test subtraction: 2 - 1 = 1.
    var diff: M.MontFieldElem = undefined;
    M.sub(&diff, sum, one_elem);
    std.debug.print("2 - 1 = {}\n", .{ M.toNormal(diff) });
    try std.testing.expect(equals(M, diff, one_elem));

    // Test negation: -1 should equal modulus - 1 (in normal representation).
    var neg: M.MontFieldElem = undefined;
    M.neg(&neg, one_elem);
    std.debug.print("Negation of 1 (normal): {}\n", .{ M.toNormal(neg) });
    try std.testing.expect(M.toNormal(neg) == (Modulus - 1));

    // Test multiplication: 2 * 3 = 6.
    const two = M.fromU32(2);
    const three = M.fromU32(3);
    var six: M.MontFieldElem = undefined;
    M.mul(&six, two, three);
    std.debug.print("2 * 3 = {}\n", .{ M.toNormal(six) });
    try std.testing.expect(M.toNormal(six) == 6);

    // Test exponentiation: 2^10 = 1024.
    var exp: M.MontFieldElem = undefined;
    M.pow(&exp, two, 10);
    std.debug.print("2^10 = {}\n", .{ M.toNormal(exp) });
    try std.testing.expect(M.toNormal(exp) == 1024);

    // Test inverse: 3 * inverse(3) should equal 1.
    var inv: M.MontFieldElem = undefined;
    M.inverse(&inv, three);
    var check: M.MontFieldElem = undefined;
    M.mul(&check, three, inv);
    std.debug.print("3 * inverse(3) = {}\n", .{ M.toNormal(check) });
    try std.testing.expect(M.toNormal(check) == 1);
}

test "MontgomeryField31 conversion to/from bytes" {
    const Modulus = 0x7fffffff;
    const M = MontgomeryField31(Modulus);

    var b: [4]u8 = undefined;
    const five = M.fromU32(5);
    M.toBytes(five, &b);
    std.debug.print("Bytes for 5: {x}\n", .{ b });
    const five_from_bytes = M.fromBytes(b[0..]);
    std.debug.print("5 recovered from bytes: {}\n", .{ M.toNormal(five_from_bytes) });
    try std.testing.expect(M.toNormal(five_from_bytes) == 5);
}

test "MontgomeryField31 conversion to/from bytes 2" {
    const Modulus = 0x7fffffff;
    const M = MontgomeryField31(Modulus);

    var b: [4]u8 = undefined;
    const value = M.fromU32(214748360);
    M.toBytes(value, &b);
    std.debug.print("Bytes : {x}\n", .{ b });
    const five_from_bytes = M.fromBytes(b[0..]);
    std.debug.print("recovered from bytes: {}\n", .{ M.toNormal(five_from_bytes) });
    try std.testing.expect(M.toNormal(five_from_bytes) == 214748360);
}