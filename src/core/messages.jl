# Error messages, wrapped for reading.
#
# The package's errors say what went wrong and what to do instead, in sentences that run to
# a few hundred characters. On one line they run off a terminal, and off the code blocks of
# the manual, which shows errors as the REPL prints them. `wrap_message` breaks long lines at
# spaces, and never inside code: a call or a backtick span stays on one line, so that it can
# be read and copied.

"""
    wrap_message(s; width = 92, first = width - 7) -> String

`s` with each line longer than `width` columns broken at spaces; the first line is held to
`first`, which leaves room for the REPL's `ERROR: `. Continuation lines keep the line's
indentation, two more for an indented line. A space inside backticks, inside brackets, or
inside the parentheses of a call such as `Interval(0, 1)` is not a break point.
"""
function wrap_message(s::AbstractString; width::Int = 92, first::Int = width - 7)
    out = String[]
    for (k, line) in enumerate(split(s, '\n'))
        limit = k == 1 ? first : width
        if textwidth(line) <= limit
            push!(out, line)
            continue
        end
        indent = match(r"^ *", line).match
        cont = isempty(indent) ? "" : indent * "  "
        cur = ""
        for w in _message_words(line)
            if isempty(cur)
                cur = indent * w
            elseif textwidth(cur) + 1 + textwidth(w) <= limit
                cur *= " " * w
            else
                push!(out, cur)
                cur, limit = cont * w, width
            end
        end
        push!(out, cur)
    end
    return join(out, '\n')
end

# The words of a line: split at spaces, except inside backticks, brackets, or the
# parentheses of a call (an opening parenthesis right after a name).
function _message_words(line::AbstractString)
    words, cur = String[], Char[]
    depth, ticks, prev = 0, false, ' '
    for c in lstrip(line, ' ')
        if c == '`'
            ticks = !ticks
        elseif !ticks && (c in ('[', '{') || (c == '(' && (depth > 0 || !isspace(prev))))
            depth += 1
        elseif !ticks && depth > 0 && c in (')', ']', '}')
            depth -= 1
        end
        if c == ' ' && !ticks && depth == 0
            isempty(cur) || push!(words, String(cur))
            empty!(cur)
        else
            push!(cur, c)
        end
        prev = c
    end
    isempty(cur) || push!(words, String(cur))
    return words
end
