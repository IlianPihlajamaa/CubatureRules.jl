# cite(rule) → BibTeX, APA or plain text (PLAN §5 item 1).

"""
    cite(rule; style = :bibtex)
    cite(citation; style = :bibtex)

The references for a rule, from its provenance: `style = :bibtex`, `:apa` or `:plain`. A
rule that records none (a user's rule passed through `transform`, say) says so rather than
returning an empty string.
"""
function cite(r::QuadratureRule; style::Symbol = :bibtex)
    isempty(r.provenance.citations) || return join((cite(c; style) for c in r.provenance.citations), "\n\n")
    style in (:bibtex, :apa, :plain) || throw(ArgumentError("unknown citation style :$style (use :bibtex, :apa or :plain)"))
    return (style === :bibtex ? "% " : "") * "no references are recorded for this rule ($(family(r)))"
end

function cite(c::Citation; style::Symbol = :bibtex)
    style === :bibtex && return _bibtex(c)
    style === :apa && return _apa(c)
    style === :plain && return _plain(c)
    throw(ArgumentError("unknown citation style :$style (use :bibtex, :apa or :plain)"))
end

function _bibtex(c::Citation)
    kind = c.kind === :software ? "misc" : c.kind === :auto ? (isempty(c.volume) && isempty(c.pages) ? "book" : "article") :
           string(c.kind)
    fields = ["author = {" * join(_bibescape.(c.authors), " and ") * "}", "title = {" * _bibescape(c.title) * "}"]
    if kind == "article"
        push!(fields, "journal = {" * _bibescape(c.journal) * "}")
    elseif kind == "misc"
        isempty(c.journal) || push!(fields, "howpublished = {" * _bibescape(c.journal) * "}")
    else
        isempty(c.journal) || push!(fields, "publisher = {" * _bibescape(c.journal) * "}")
    end
    push!(fields, "year = {" * string(c.year) * "}")
    isempty(c.volume) || push!(fields, "volume = {" * c.volume * "}")
    isempty(c.pages) || push!(fields, "pages = {" * c.pages * "}")
    isempty(c.doi) || push!(fields, "doi = {" * c.doi * "}")
    isempty(c.url) || push!(fields, "url = {" * c.url * "}")
    return "@" * kind * "{" * c.key * ",\n  " * join(fields, ",\n  ") * "\n}"
end

# the characters BibTeX would otherwise read as markup
_bibescape(s::AbstractString) = replace(s, r"([&%$#_])" => s"\\\1")

function _apa_name(n)
    parts = split(n)
    length(parts) == 1 && return n
    initials = join((string(first(p), ".") for p in parts[1:(end - 1)]), " ")
    return parts[end] * ", " * initials
end

function _apa(c::Citation)
    names = _apa_name.(c.authors)
    auth = length(names) == 1 ? names[1] :
           join(names[1:(end - 1)], ", ") * ", & " * names[end]
    s = auth * " (" * string(c.year) * "). " * c.title * ". " * c.journal
    isempty(c.volume) || (s *= ", " * c.volume)
    isempty(c.pages) || (s *= ", " * replace(c.pages, "--" => "–"))
    s *= "."
    isempty(c.doi) || (s *= " https://doi.org/" * c.doi)
    isempty(c.doi) && !isempty(c.url) && (s *= " " * c.url)
    return s
end

function _plain(c::Citation)
    s = join(c.authors, ", ") * ", \"" * c.title * "\", " * c.journal
    isempty(c.volume) || (s *= " " * c.volume)
    s *= " (" * string(c.year) * ")"
    isempty(c.pages) || (s *= " " * replace(c.pages, "--" => "–"))
    isempty(c.doi) || (s *= ", doi:" * c.doi)
    isempty(c.doi) && !isempty(c.url) && (s *= ", " * c.url)
    return s
end
