# Sanitize Excel worksheet names

Ensures worksheet names are compatible with Excel and `openxlsx`. Names
are:

- Truncated to `max_len` characters (default 31, Excel's limit).

- Cleaned of invalid characters (`: \ / ? * [ ]`).

- Made unique by appending numeric suffixes (e.g. `"name"`, `"name_1"`).

If any names are modified, a warning is emitted listing the original and
sanitized names.

## Usage

``` r
sanitize_sheet_names(names, max_len = 31L)
```

## Arguments

- names:

  Character vector of proposed worksheet names.

- max_len:

  Integer maximum length for each name (default 31).

## Value

A character vector of sanitized, unique worksheet names.
