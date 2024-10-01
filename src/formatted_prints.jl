"""
    @File: formatted_prints.jl
    @Author: Jindriska Deckerova, Petr Vana, Jan Faigl
"""

using Printf

function print_info(msg::String)
    printstyled(stdout, "Info: ", color=:blue, bold=true)
    printstyled(stdout, "$(msg)\n")
end

function print_update(msg::String)
    printstyled(stdout, "Update: ", color=:magenta, bold=true)
    printstyled(stdout, "$(msg)\n")
end

function print_plain_message(msg::String, color)
    printstyled(stdout, "$(msg)\n", color=color, bold=false)
end

function print_success_message(msg::String)
    printstyled(stdout, "SUCCESS: ", color=:green, bold=true)
    printstyled(stdout, "$(msg)\n", color=:green, bold=false)
end

function print_fail_message(msg::String)
    printstyled(stdout, "FAIL: ", color=:red, bold=true)
    printstyled(stdout, "$(msg)\n")
end

function print_warn_message(msg::String)
    printstyled(stdout, "WARN: ", color=:yellow, bold=true)
    printstyled(stdout, "$(msg)\n")
end