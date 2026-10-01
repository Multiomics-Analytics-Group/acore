# acore.types package

Collect common types of pandas DataFrames used in the package.

Documentation of DataFrame Models API:
[https://pandera.readthedocs.io/en/stable/dataframe_models.html](https://pandera.readthedocs.io/en/stable/dataframe_models.html)

### check_numeric_dataframe(df: DataFrame) → DataFrame

Check if the DataFrame contains only numeric data.
returns the DataFrame again if it is valid (allowing chaining).

### select_numeric_columns(df: DataFrame) → DataFrame

Select only numeric columns from the DataFrame.

### build_schema_all_floats(df: DataFrame) → [DataFrameSchema](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.api.pandas.container.DataFrameSchema.html#pandera.api.pandas.container.DataFrameSchema)

Build a schema that checks if all columns are float, potentially
containing NaN values.

## Submodules

## acore.types.differential_analysis module

### *class* AnovaSchema(\*args, \*\*kwargs)

Bases: [`DataFrameModel`](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.api.pandas.model.DataFrameModel.html#pandera.api.pandas.model.DataFrameModel)

Schema for the enrichment analysis results DataFrame.

#### group1 *: str* *= 'group1'*

#### group2 *: str* *= 'group2'*

#### mean_group1 *: float* *= 'mean(group1)'*

#### std_group1 *: float* *= 'std(group1)'*

#### mean_group2 *: float* *= 'mean(group2)'*

#### std_group2 *: float* *= 'std(group2)'*

#### t_statistics *: float* *= 'T-Statistics'*

#### pvalue *: float* *= 'pvalue'*

#### log2FC *: float* *= 'log2FC'*

#### FC *: float* *= 'FC'*

#### padj *: float* *= 'padj'*

#### correction *: str* *= 'correction'*

#### rejected *: bool* *= 'rejected'*

#### neg_log10_p_value *: float* *= '-log10 pvalue'*

#### Method *: str* *= 'Method'*

#### *class* Config

Bases: [`BaseConfig`](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.api.pandas.model_config.BaseConfig.html#pandera.api.pandas.model_config.BaseConfig)

#### add_missing_columns *: bool* *= False*

add columns to dataframe if they are missing

#### coerce *: bool* *= False*

coerce types of all schema components

#### description *: str | None* *= None*

arbitrary textual description

#### drop_invalid_rows *: bool* *= False*

drop invalid rows on validation

#### dtype *: str | type | [DataType](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.dtypes.DataType.html#pandera.dtypes.DataType) | ExtensionDtype | [dtype](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype) | None* *= None*

datatype of the dataframe. This overrides the data types specified in
any of the fields.

#### from_format *: Literal[Formats.csv] | Literal[Formats.dict] | Literal[Formats.json] | Literal[Formats.feather] | Literal[Formats.parquet] | Literal[Formats.pickle] | Literal[Formats.json_normalize] | Callable | None* *= None*

data format before validation. This option only applies to
schemas used in the context of the pandera type constructor
`pa.typing.DataFrame[Schema](data)`. If None, assumes a data structure
compatible with the `pandas.DataFrame` constructor.

#### from_format_kwargs *: dict[str, Any] | None* *= None*

a dictionary keyword arguments to pass into the reader function that
converts the object of type `from_format` to a pandera-validate-able
data structure. The reader function is implemented in the pandera.typing
generic types via the `from_format` and `to_format` methods.

#### metadata *: dict | None* *= None*

a dictionary object to store key-value data at schema level

#### multiindex_coerce *: bool* *= False*

coerce types of all MultiIndex components

#### multiindex_name *: str | None* *= None*

name of multiindex

#### multiindex_ordered *: bool* *= True*

validate MultiIndex in order

#### multiindex_strict *: bool | Literal['filter']* *= False*

make sure all specified columns are in validated MultiIndex -
if `"filter"`, removes indexes not specified in the schema

#### multiindex_unique *= None*

make sure the MultiIndex is unique along the list of columns

#### name *: str | None* *= 'AnovaSchema'*

name of schema

#### on_missing_columns *: str | None* *= None*

action to take when an optional column is missing from the dataframe.
`None` (default) passes silently, `"warn"` emits a
`SchemaWarning`. Overridable per field via
`Field(on_missing=...)`.

#### ordered *: bool* *= False*

validate columns order

#### strict *: bool | Literal['filter']* *= False*

make sure all specified columns are in the validated dataframe -
if `"filter"`, removes columns not specified in the schema

#### title *: str | None* *= None*

human-readable label for schema

#### to_format *: Literal[Formats.csv] | Literal[Formats.dict] | Literal[Formats.json] | Literal[Formats.feather] | Literal[Formats.parquet] | Literal[Formats.pickle] | Literal[Formats.json_normalize] | Callable | None* *= None*

data format to serialize into after validation. This option only applies
to  schemas used in the context of the pandera type constructor
`pa.typing.DataFrame[Schema](data)`. If None, returns a dataframe.

#### to_format_buffer *: str | Callable | None* *= None*

Buffer to be provided when to_format is a custom callable. See docs for
example of how to implement an example of a to format function.

#### to_format_kwargs *: dict[str, Any] | None* *= None*

a dictionary keyword arguments to pass into the writer function that
converts the pandera-validate-able object to type `to_format`.
The writer function is implemented in the pandera.typing
generic types via the `from_format` and `to_format` methods.

#### unique *: str | list[str] | None* *= None*

make sure certain column combinations are unique

#### unique_column_names *: bool* *= False*

make sure dataframe column names are unique

#### *classmethod* build_schema_(\*\*kwargs) → [DataFrameSchema](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.api.pandas.container.DataFrameSchema.html#pandera.api.pandas.container.DataFrameSchema)

#### *classmethod* empty(\*\_args) → [DataFrame](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.typing.DataFrame.html#pandera.typing.DataFrame)[Self]

Create an empty DataFrame with the schema of this model.

#### *classmethod* example(\*\*kwargs) → DataFrameBase[Self]

Generate an example of this data model specification.

#### *classmethod* from_json(source)

Load a schema from JSON.

* **Parameters:**
  **source** – str, Path, or file stream with JSON content.
* **Returns:**
  the backend-specific schema object.

#### *classmethod* from_yaml(yaml_schema)

Load a schema from YAML.

* **Parameters:**
  **yaml_schema** – str, Path, or file stream with YAML content.
* **Returns:**
  the backend-specific schema object.

#### *classmethod* get_metadata() → dict | None

Provide metadata for columns and schema level

#### *classmethod* pydantic_validate(schema_model: Any) → [DataFrameModel](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.api.dataframe.model.DataFrameModel.html#pandera.api.dataframe.model.DataFrameModel)

Verify that the input is a compatible dataframe model.

#### *classmethod* strategy(\*\*kwargs)

Create a data synthesis strategy.

#### *classmethod* to_json(target: PathLike | None = None, \*\*kwargs)

Convert this model’s schema to JSON.

#### *classmethod* to_json_schema()

Serialize schema metadata into json-schema format.

* **Parameters:**
  **dataframe_schema** – schema to write to json-schema format.

#### NOTE
This function is currently does not fully specify a pandera schema,
and is primarily used internally to render OpenAPI docs via the
FastAPI integration.

#### *classmethod* to_schema() → TSchema

Create `DataFrameSchema` from the `DataFrameModel`.

#### *classmethod* to_yaml(stream: PathLike | None = None)

Convert this model’s schema to YAML.

#### *classmethod* validate(check_obj: DataFrame, head: int | None = None, tail: int | None = None, sample: int | None = None, random_state: int | None = None, lazy: bool = False, inplace: bool = False) → [DataFrame](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.typing.DataFrame.html#pandera.typing.DataFrame)[Self]

Validate a DataFrame based on the schema specification.

* **Parameters:**
  * **check_obj** (*pd.DataFrame*) – the dataframe to be validated.
  * **head** – validate the first n rows. Rows overlapping with tail or
    sample are de-duplicated.
  * **tail** – validate the last n rows. Rows overlapping with head or
    sample are de-duplicated.
  * **sample** – validate a random sample of n rows. Rows overlapping
    with head or tail are de-duplicated.
  * **random_state** – random seed for the `sample` argument.
  * **lazy** – if True, lazily evaluates dataframe against all validation
    checks and raises a `SchemaErrors`. Otherwise, raise
    `SchemaError` as soon as one occurs.
  * **inplace** – if True, applies coercion to the object of validation,
    otherwise creates a copy of the data.
* **Returns:**
  validated `DataFrame`
* **Raises:**
  **SchemaError** – when `DataFrame` violates built-in or custom
  checks.

### *class* AnovaSchemaMultiGroup(\*args, \*\*kwargs)

Bases: [`AnovaSchema`](#acore.types.differential_analysis.AnovaSchema)

Schema for more than two groups

#### t_statistics *: float* *= 'posthoc T-Statistics'*

#### posthoc_pvalue *: float* *= 'posthoc pvalue'*

#### f_statistics *: float* *= 'F-statistics'*

#### posthoc_padj *: float* *= 'posthoc padj'*

#### posthoc_paired *: bool* *= 'posthoc Paired'*

#### posthoc_parametric *: bool* *= 'posthoc Parametric'*

#### posthoc_dof *: float* *= 'posthoc dof'*

#### posthoc_tail *: str* *= 'posthoc tail'*

#### posthoc_BF10 *: float* *= 'posthoc BF10'*

#### posthoc_effsize *: float* *= 'posthoc effsize'*

#### efftype *: str* *= 'efftype'*

#### *class* Config

Bases: [`Config`](#acore.types.exploratory_analysis.TwoLoadingsSchema.Config)

#### add_missing_columns *: bool* *= False*

add columns to dataframe if they are missing

#### coerce *: bool* *= False*

coerce types of all schema components

#### description *: str | None* *= None*

arbitrary textual description

#### drop_invalid_rows *: bool* *= False*

drop invalid rows on validation

#### dtype *: str | type | [DataType](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.dtypes.DataType.html#pandera.dtypes.DataType) | ExtensionDtype | [dtype](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype) | None* *= None*

datatype of the dataframe. This overrides the data types specified in
any of the fields.

#### from_format *: Literal[Formats.csv] | Literal[Formats.dict] | Literal[Formats.json] | Literal[Formats.feather] | Literal[Formats.parquet] | Literal[Formats.pickle] | Literal[Formats.json_normalize] | Callable | None* *= None*

data format before validation. This option only applies to
schemas used in the context of the pandera type constructor
`pa.typing.DataFrame[Schema](data)`. If None, assumes a data structure
compatible with the `pandas.DataFrame` constructor.

#### from_format_kwargs *: dict[str, Any] | None* *= None*

a dictionary keyword arguments to pass into the reader function that
converts the object of type `from_format` to a pandera-validate-able
data structure. The reader function is implemented in the pandera.typing
generic types via the `from_format` and `to_format` methods.

#### metadata *: dict | None* *= None*

a dictionary object to store key-value data at schema level

#### multiindex_coerce *: bool* *= False*

coerce types of all MultiIndex components

#### multiindex_name *: str | None* *= None*

name of multiindex

#### multiindex_ordered *: bool* *= True*

validate MultiIndex in order

#### multiindex_strict *: bool | Literal['filter']* *= False*

make sure all specified columns are in validated MultiIndex -
if `"filter"`, removes indexes not specified in the schema

#### multiindex_unique *= None*

make sure the MultiIndex is unique along the list of columns

#### name *: str | None* *= 'AnovaSchemaMultiGroup'*

name of schema

#### on_missing_columns *: str | None* *= None*

action to take when an optional column is missing from the dataframe.
`None` (default) passes silently, `"warn"` emits a
`SchemaWarning`. Overridable per field via
`Field(on_missing=...)`.

#### ordered *: bool* *= False*

validate columns order

#### strict *: bool | Literal['filter']* *= False*

make sure all specified columns are in the validated dataframe -
if `"filter"`, removes columns not specified in the schema

#### title *: str | None* *= None*

human-readable label for schema

#### to_format *: Literal[Formats.csv] | Literal[Formats.dict] | Literal[Formats.json] | Literal[Formats.feather] | Literal[Formats.parquet] | Literal[Formats.pickle] | Literal[Formats.json_normalize] | Callable | None* *= None*

data format to serialize into after validation. This option only applies
to  schemas used in the context of the pandera type constructor
`pa.typing.DataFrame[Schema](data)`. If None, returns a dataframe.

#### to_format_buffer *: str | Callable | None* *= None*

Buffer to be provided when to_format is a custom callable. See docs for
example of how to implement an example of a to format function.

#### to_format_kwargs *: dict[str, Any] | None* *= None*

a dictionary keyword arguments to pass into the writer function that
converts the pandera-validate-able object to type `to_format`.
The writer function is implemented in the pandera.typing
generic types via the `from_format` and `to_format` methods.

#### unique *: str | list[str] | None* *= None*

make sure certain column combinations are unique

#### unique_column_names *: bool* *= False*

make sure dataframe column names are unique

#### FC *: float* *= 'FC'*

#### Method *: str* *= 'Method'*

#### *classmethod* build_schema_(\*\*kwargs) → [DataFrameSchema](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.api.pandas.container.DataFrameSchema.html#pandera.api.pandas.container.DataFrameSchema)

#### correction *: str* *= 'correction'*

#### *classmethod* empty(\*\_args) → [DataFrame](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.typing.DataFrame.html#pandera.typing.DataFrame)[Self]

Create an empty DataFrame with the schema of this model.

#### *classmethod* example(\*\*kwargs) → DataFrameBase[Self]

Generate an example of this data model specification.

#### *classmethod* from_json(source)

Load a schema from JSON.

* **Parameters:**
  **source** – str, Path, or file stream with JSON content.
* **Returns:**
  the backend-specific schema object.

#### *classmethod* from_yaml(yaml_schema)

Load a schema from YAML.

* **Parameters:**
  **yaml_schema** – str, Path, or file stream with YAML content.
* **Returns:**
  the backend-specific schema object.

#### *classmethod* get_metadata() → dict | None

Provide metadata for columns and schema level

#### group1 *: str* *= 'group1'*

#### group2 *: str* *= 'group2'*

#### log2FC *: float* *= 'log2FC'*

#### mean_group1 *: float* *= 'mean(group1)'*

#### mean_group2 *: float* *= 'mean(group2)'*

#### neg_log10_p_value *: float* *= '-log10 pvalue'*

#### padj *: float* *= 'padj'*

#### pvalue *: float* *= 'pvalue'*

#### *classmethod* pydantic_validate(schema_model: Any) → [DataFrameModel](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.api.dataframe.model.DataFrameModel.html#pandera.api.dataframe.model.DataFrameModel)

Verify that the input is a compatible dataframe model.

#### rejected *: bool* *= 'rejected'*

#### std_group1 *: float* *= 'std(group1)'*

#### std_group2 *: float* *= 'std(group2)'*

#### *classmethod* strategy(\*\*kwargs)

Create a data synthesis strategy.

#### *classmethod* to_json(target: PathLike | None = None, \*\*kwargs)

Convert this model’s schema to JSON.

#### *classmethod* to_json_schema()

Serialize schema metadata into json-schema format.

* **Parameters:**
  **dataframe_schema** – schema to write to json-schema format.

#### NOTE
This function is currently does not fully specify a pandera schema,
and is primarily used internally to render OpenAPI docs via the
FastAPI integration.

#### *classmethod* to_schema() → TSchema

Create `DataFrameSchema` from the `DataFrameModel`.

#### *classmethod* to_yaml(stream: PathLike | None = None)

Convert this model’s schema to YAML.

#### *classmethod* validate(check_obj: DataFrame, head: int | None = None, tail: int | None = None, sample: int | None = None, random_state: int | None = None, lazy: bool = False, inplace: bool = False) → [DataFrame](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.typing.DataFrame.html#pandera.typing.DataFrame)[Self]

Validate a DataFrame based on the schema specification.

* **Parameters:**
  * **check_obj** (*pd.DataFrame*) – the dataframe to be validated.
  * **head** – validate the first n rows. Rows overlapping with tail or
    sample are de-duplicated.
  * **tail** – validate the last n rows. Rows overlapping with head or
    sample are de-duplicated.
  * **sample** – validate a random sample of n rows. Rows overlapping
    with head or tail are de-duplicated.
  * **random_state** – random seed for the `sample` argument.
  * **lazy** – if True, lazily evaluates dataframe against all validation
    checks and raises a `SchemaErrors`. Otherwise, raise
    `SchemaError` as soon as one occurs.
  * **inplace** – if True, applies coercion to the object of validation,
    otherwise creates a copy of the data.
* **Returns:**
  validated `DataFrame`
* **Raises:**
  **SchemaError** – when `DataFrame` violates built-in or custom
  checks.

### *class* AncovaSchema(\*args, \*\*kwargs)

Bases: [`DataFrameModel`](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.api.pandas.model.DataFrameModel.html#pandera.api.pandas.model.DataFrameModel)

Schema for the enrichment analysis results DataFrame.

#### t_statistics *: float* *= 'posthoc T-Statistics'*

#### posthoc_pvalue *: float* *= 'posthoc pvalue'*

#### coef *: float* *= 'coef'*

#### std_err *: float* *= 'std err'*

#### conf_int_low *: float* *= 'Conf. Int. Low'*

#### conf_int_upp *: float* *= 'Conf. Int. Upp.'*

#### f_statistics *: float* *= 'F-statistics'*

#### posthoc_padj *: float* *= 'posthoc padj'*

#### *class* Config

Bases: [`BaseConfig`](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.api.pandas.model_config.BaseConfig.html#pandera.api.pandas.model_config.BaseConfig)

#### add_missing_columns *: bool* *= False*

add columns to dataframe if they are missing

#### coerce *: bool* *= False*

coerce types of all schema components

#### description *: str | None* *= None*

arbitrary textual description

#### drop_invalid_rows *: bool* *= False*

drop invalid rows on validation

#### dtype *: str | type | [DataType](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.dtypes.DataType.html#pandera.dtypes.DataType) | ExtensionDtype | [dtype](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype) | None* *= None*

datatype of the dataframe. This overrides the data types specified in
any of the fields.

#### from_format *: Literal[Formats.csv] | Literal[Formats.dict] | Literal[Formats.json] | Literal[Formats.feather] | Literal[Formats.parquet] | Literal[Formats.pickle] | Literal[Formats.json_normalize] | Callable | None* *= None*

data format before validation. This option only applies to
schemas used in the context of the pandera type constructor
`pa.typing.DataFrame[Schema](data)`. If None, assumes a data structure
compatible with the `pandas.DataFrame` constructor.

#### from_format_kwargs *: dict[str, Any] | None* *= None*

a dictionary keyword arguments to pass into the reader function that
converts the object of type `from_format` to a pandera-validate-able
data structure. The reader function is implemented in the pandera.typing
generic types via the `from_format` and `to_format` methods.

#### metadata *: dict | None* *= None*

a dictionary object to store key-value data at schema level

#### multiindex_coerce *: bool* *= False*

coerce types of all MultiIndex components

#### multiindex_name *: str | None* *= None*

name of multiindex

#### multiindex_ordered *: bool* *= True*

validate MultiIndex in order

#### multiindex_strict *: bool | Literal['filter']* *= False*

make sure all specified columns are in validated MultiIndex -
if `"filter"`, removes indexes not specified in the schema

#### multiindex_unique *= None*

make sure the MultiIndex is unique along the list of columns

#### name *: str | None* *= 'AncovaSchema'*

name of schema

#### on_missing_columns *: str | None* *= None*

action to take when an optional column is missing from the dataframe.
`None` (default) passes silently, `"warn"` emits a
`SchemaWarning`. Overridable per field via
`Field(on_missing=...)`.

#### ordered *: bool* *= False*

validate columns order

#### strict *: bool | Literal['filter']* *= False*

make sure all specified columns are in the validated dataframe -
if `"filter"`, removes columns not specified in the schema

#### title *: str | None* *= None*

human-readable label for schema

#### to_format *: Literal[Formats.csv] | Literal[Formats.dict] | Literal[Formats.json] | Literal[Formats.feather] | Literal[Formats.parquet] | Literal[Formats.pickle] | Literal[Formats.json_normalize] | Callable | None* *= None*

data format to serialize into after validation. This option only applies
to  schemas used in the context of the pandera type constructor
`pa.typing.DataFrame[Schema](data)`. If None, returns a dataframe.

#### to_format_buffer *: str | Callable | None* *= None*

Buffer to be provided when to_format is a custom callable. See docs for
example of how to implement an example of a to format function.

#### to_format_kwargs *: dict[str, Any] | None* *= None*

a dictionary keyword arguments to pass into the writer function that
converts the pandera-validate-able object to type `to_format`.
The writer function is implemented in the pandera.typing
generic types via the `from_format` and `to_format` methods.

#### unique *: str | list[str] | None* *= None*

make sure certain column combinations are unique

#### unique_column_names *: bool* *= False*

make sure dataframe column names are unique

#### *classmethod* build_schema_(\*\*kwargs) → [DataFrameSchema](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.api.pandas.container.DataFrameSchema.html#pandera.api.pandas.container.DataFrameSchema)

#### *classmethod* empty(\*\_args) → [DataFrame](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.typing.DataFrame.html#pandera.typing.DataFrame)[Self]

Create an empty DataFrame with the schema of this model.

#### *classmethod* example(\*\*kwargs) → DataFrameBase[Self]

Generate an example of this data model specification.

#### *classmethod* from_json(source)

Load a schema from JSON.

* **Parameters:**
  **source** – str, Path, or file stream with JSON content.
* **Returns:**
  the backend-specific schema object.

#### *classmethod* from_yaml(yaml_schema)

Load a schema from YAML.

* **Parameters:**
  **yaml_schema** – str, Path, or file stream with YAML content.
* **Returns:**
  the backend-specific schema object.

#### *classmethod* get_metadata() → dict | None

Provide metadata for columns and schema level

#### *classmethod* pydantic_validate(schema_model: Any) → [DataFrameModel](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.api.dataframe.model.DataFrameModel.html#pandera.api.dataframe.model.DataFrameModel)

Verify that the input is a compatible dataframe model.

#### *classmethod* strategy(\*\*kwargs)

Create a data synthesis strategy.

#### *classmethod* to_json(target: PathLike | None = None, \*\*kwargs)

Convert this model’s schema to JSON.

#### *classmethod* to_json_schema()

Serialize schema metadata into json-schema format.

* **Parameters:**
  **dataframe_schema** – schema to write to json-schema format.

#### NOTE
This function is currently does not fully specify a pandera schema,
and is primarily used internally to render OpenAPI docs via the
FastAPI integration.

#### *classmethod* to_schema() → TSchema

Create `DataFrameSchema` from the `DataFrameModel`.

#### *classmethod* to_yaml(stream: PathLike | None = None)

Convert this model’s schema to YAML.

#### *classmethod* validate(check_obj: DataFrame, head: int | None = None, tail: int | None = None, sample: int | None = None, random_state: int | None = None, lazy: bool = False, inplace: bool = False) → [DataFrame](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.typing.DataFrame.html#pandera.typing.DataFrame)[Self]

Validate a DataFrame based on the schema specification.

* **Parameters:**
  * **check_obj** (*pd.DataFrame*) – the dataframe to be validated.
  * **head** – validate the first n rows. Rows overlapping with tail or
    sample are de-duplicated.
  * **tail** – validate the last n rows. Rows overlapping with head or
    sample are de-duplicated.
  * **sample** – validate a random sample of n rows. Rows overlapping
    with head or tail are de-duplicated.
  * **random_state** – random seed for the `sample` argument.
  * **lazy** – if True, lazily evaluates dataframe against all validation
    checks and raises a `SchemaErrors`. Otherwise, raise
    `SchemaError` as soon as one occurs.
  * **inplace** – if True, applies coercion to the object of validation,
    otherwise creates a copy of the data.
* **Returns:**
  validated `DataFrame`
* **Raises:**
  **SchemaError** – when `DataFrame` violates built-in or custom
  checks.

## acore.types.enrichment_analysis module

### *class* EnrichmentAnalysisSchema(\*args, \*\*kwargs)

Bases: [`DataFrameModel`](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.api.pandas.model.DataFrameModel.html#pandera.api.pandas.model.DataFrameModel)

Schema for the enrichment analysis results DataFrame.

#### terms *: str* *= 'terms'*

#### identifiers *: str* *= 'identifiers'*

#### foreground *: int* *= 'foreground'*

#### background *: int* *= 'background'*

#### foreground_pop *: int* *= 'foreground_pop'*

#### background_pop *: int* *= 'background_pop'*

#### pvalue *: float* *= 'pvalue'*

#### padj *: float* *= 'padj'*

#### rejected *: bool* *= 'rejected'*

#### direction *: str* *= 'direction'*

#### comparison *: str* *= 'comparison'*

#### *class* Config

Bases: [`BaseConfig`](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.api.pandas.model_config.BaseConfig.html#pandera.api.pandas.model_config.BaseConfig)

#### add_missing_columns *: bool* *= False*

add columns to dataframe if they are missing

#### coerce *: bool* *= False*

coerce types of all schema components

#### description *: str | None* *= None*

arbitrary textual description

#### drop_invalid_rows *: bool* *= False*

drop invalid rows on validation

#### dtype *: str | type | [DataType](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.dtypes.DataType.html#pandera.dtypes.DataType) | ExtensionDtype | [dtype](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype) | None* *= None*

datatype of the dataframe. This overrides the data types specified in
any of the fields.

#### from_format *: Literal[Formats.csv] | Literal[Formats.dict] | Literal[Formats.json] | Literal[Formats.feather] | Literal[Formats.parquet] | Literal[Formats.pickle] | Literal[Formats.json_normalize] | Callable | None* *= None*

data format before validation. This option only applies to
schemas used in the context of the pandera type constructor
`pa.typing.DataFrame[Schema](data)`. If None, assumes a data structure
compatible with the `pandas.DataFrame` constructor.

#### from_format_kwargs *: dict[str, Any] | None* *= None*

a dictionary keyword arguments to pass into the reader function that
converts the object of type `from_format` to a pandera-validate-able
data structure. The reader function is implemented in the pandera.typing
generic types via the `from_format` and `to_format` methods.

#### metadata *: dict | None* *= None*

a dictionary object to store key-value data at schema level

#### multiindex_coerce *: bool* *= False*

coerce types of all MultiIndex components

#### multiindex_name *: str | None* *= None*

name of multiindex

#### multiindex_ordered *: bool* *= True*

validate MultiIndex in order

#### multiindex_strict *: bool | Literal['filter']* *= False*

make sure all specified columns are in validated MultiIndex -
if `"filter"`, removes indexes not specified in the schema

#### multiindex_unique *= None*

make sure the MultiIndex is unique along the list of columns

#### name *: str | None* *= 'EnrichmentAnalysisSchema'*

name of schema

#### on_missing_columns *: str | None* *= None*

action to take when an optional column is missing from the dataframe.
`None` (default) passes silently, `"warn"` emits a
`SchemaWarning`. Overridable per field via
`Field(on_missing=...)`.

#### ordered *: bool* *= False*

validate columns order

#### strict *: bool | Literal['filter']* *= False*

make sure all specified columns are in the validated dataframe -
if `"filter"`, removes columns not specified in the schema

#### title *: str | None* *= None*

human-readable label for schema

#### to_format *: Literal[Formats.csv] | Literal[Formats.dict] | Literal[Formats.json] | Literal[Formats.feather] | Literal[Formats.parquet] | Literal[Formats.pickle] | Literal[Formats.json_normalize] | Callable | None* *= None*

data format to serialize into after validation. This option only applies
to  schemas used in the context of the pandera type constructor
`pa.typing.DataFrame[Schema](data)`. If None, returns a dataframe.

#### to_format_buffer *: str | Callable | None* *= None*

Buffer to be provided when to_format is a custom callable. See docs for
example of how to implement an example of a to format function.

#### to_format_kwargs *: dict[str, Any] | None* *= None*

a dictionary keyword arguments to pass into the writer function that
converts the pandera-validate-able object to type `to_format`.
The writer function is implemented in the pandera.typing
generic types via the `from_format` and `to_format` methods.

#### unique *: str | list[str] | None* *= None*

make sure certain column combinations are unique

#### unique_column_names *: bool* *= False*

make sure dataframe column names are unique

#### *classmethod* build_schema_(\*\*kwargs) → [DataFrameSchema](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.api.pandas.container.DataFrameSchema.html#pandera.api.pandas.container.DataFrameSchema)

#### *classmethod* empty(\*\_args) → [DataFrame](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.typing.DataFrame.html#pandera.typing.DataFrame)[Self]

Create an empty DataFrame with the schema of this model.

#### *classmethod* example(\*\*kwargs) → DataFrameBase[Self]

Generate an example of this data model specification.

#### *classmethod* from_json(source)

Load a schema from JSON.

* **Parameters:**
  **source** – str, Path, or file stream with JSON content.
* **Returns:**
  the backend-specific schema object.

#### *classmethod* from_yaml(yaml_schema)

Load a schema from YAML.

* **Parameters:**
  **yaml_schema** – str, Path, or file stream with YAML content.
* **Returns:**
  the backend-specific schema object.

#### *classmethod* get_metadata() → dict | None

Provide metadata for columns and schema level

#### *classmethod* pydantic_validate(schema_model: Any) → [DataFrameModel](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.api.dataframe.model.DataFrameModel.html#pandera.api.dataframe.model.DataFrameModel)

Verify that the input is a compatible dataframe model.

#### *classmethod* strategy(\*\*kwargs)

Create a data synthesis strategy.

#### *classmethod* to_json(target: PathLike | None = None, \*\*kwargs)

Convert this model’s schema to JSON.

#### *classmethod* to_json_schema()

Serialize schema metadata into json-schema format.

* **Parameters:**
  **dataframe_schema** – schema to write to json-schema format.

#### NOTE
This function is currently does not fully specify a pandera schema,
and is primarily used internally to render OpenAPI docs via the
FastAPI integration.

#### *classmethod* to_schema() → TSchema

Create `DataFrameSchema` from the `DataFrameModel`.

#### *classmethod* to_yaml(stream: PathLike | None = None)

Convert this model’s schema to YAML.

#### *classmethod* validate(check_obj: DataFrame, head: int | None = None, tail: int | None = None, sample: int | None = None, random_state: int | None = None, lazy: bool = False, inplace: bool = False) → [DataFrame](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.typing.DataFrame.html#pandera.typing.DataFrame)[Self]

Validate a DataFrame based on the schema specification.

* **Parameters:**
  * **check_obj** (*pd.DataFrame*) – the dataframe to be validated.
  * **head** – validate the first n rows. Rows overlapping with tail or
    sample are de-duplicated.
  * **tail** – validate the last n rows. Rows overlapping with head or
    sample are de-duplicated.
  * **sample** – validate a random sample of n rows. Rows overlapping
    with head or tail are de-duplicated.
  * **random_state** – random seed for the `sample` argument.
  * **lazy** – if True, lazily evaluates dataframe against all validation
    checks and raises a `SchemaErrors`. Otherwise, raise
    `SchemaError` as soon as one occurs.
  * **inplace** – if True, applies coercion to the object of validation,
    otherwise creates a copy of the data.
* **Returns:**
  validated `DataFrame`
* **Raises:**
  **SchemaError** – when `DataFrame` violates built-in or custom
  checks.

## acore.types.exploratory_analysis module

### *class* TwoComponentSchema(\*args, \*\*kwargs)

Bases: [`DataFrameModel`](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.api.pandas.model.DataFrameModel.html#pandera.api.pandas.model.DataFrameModel)

Schema for the PCA components DataFrame.

#### group *: str | None* *= 'group'*

#### x *: float* *= 'x'*

#### y *: float* *= 'y'*

#### *class* Config

Bases: [`BaseConfig`](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.api.pandas.model_config.BaseConfig.html#pandera.api.pandas.model_config.BaseConfig)

#### add_missing_columns *: bool* *= False*

add columns to dataframe if they are missing

#### coerce *: bool* *= False*

coerce types of all schema components

#### description *: str | None* *= None*

arbitrary textual description

#### drop_invalid_rows *: bool* *= False*

drop invalid rows on validation

#### dtype *: str | type | [DataType](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.dtypes.DataType.html#pandera.dtypes.DataType) | ExtensionDtype | [dtype](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype) | None* *= None*

datatype of the dataframe. This overrides the data types specified in
any of the fields.

#### from_format *: Literal[Formats.csv] | Literal[Formats.dict] | Literal[Formats.json] | Literal[Formats.feather] | Literal[Formats.parquet] | Literal[Formats.pickle] | Literal[Formats.json_normalize] | Callable | None* *= None*

data format before validation. This option only applies to
schemas used in the context of the pandera type constructor
`pa.typing.DataFrame[Schema](data)`. If None, assumes a data structure
compatible with the `pandas.DataFrame` constructor.

#### from_format_kwargs *: dict[str, Any] | None* *= None*

a dictionary keyword arguments to pass into the reader function that
converts the object of type `from_format` to a pandera-validate-able
data structure. The reader function is implemented in the pandera.typing
generic types via the `from_format` and `to_format` methods.

#### metadata *: dict | None* *= None*

a dictionary object to store key-value data at schema level

#### multiindex_coerce *: bool* *= False*

coerce types of all MultiIndex components

#### multiindex_name *: str | None* *= None*

name of multiindex

#### multiindex_ordered *: bool* *= True*

validate MultiIndex in order

#### multiindex_strict *: bool | Literal['filter']* *= False*

make sure all specified columns are in validated MultiIndex -
if `"filter"`, removes indexes not specified in the schema

#### multiindex_unique *= None*

make sure the MultiIndex is unique along the list of columns

#### name *: str | None* *= 'TwoComponentSchema'*

name of schema

#### on_missing_columns *: str | None* *= None*

action to take when an optional column is missing from the dataframe.
`None` (default) passes silently, `"warn"` emits a
`SchemaWarning`. Overridable per field via
`Field(on_missing=...)`.

#### ordered *: bool* *= False*

validate columns order

#### strict *: bool | Literal['filter']* *= False*

make sure all specified columns are in the validated dataframe -
if `"filter"`, removes columns not specified in the schema

#### title *: str | None* *= None*

human-readable label for schema

#### to_format *: Literal[Formats.csv] | Literal[Formats.dict] | Literal[Formats.json] | Literal[Formats.feather] | Literal[Formats.parquet] | Literal[Formats.pickle] | Literal[Formats.json_normalize] | Callable | None* *= None*

data format to serialize into after validation. This option only applies
to  schemas used in the context of the pandera type constructor
`pa.typing.DataFrame[Schema](data)`. If None, returns a dataframe.

#### to_format_buffer *: str | Callable | None* *= None*

Buffer to be provided when to_format is a custom callable. See docs for
example of how to implement an example of a to format function.

#### to_format_kwargs *: dict[str, Any] | None* *= None*

a dictionary keyword arguments to pass into the writer function that
converts the pandera-validate-able object to type `to_format`.
The writer function is implemented in the pandera.typing
generic types via the `from_format` and `to_format` methods.

#### unique *: str | list[str] | None* *= None*

make sure certain column combinations are unique

#### unique_column_names *: bool* *= False*

make sure dataframe column names are unique

#### *classmethod* build_schema_(\*\*kwargs) → [DataFrameSchema](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.api.pandas.container.DataFrameSchema.html#pandera.api.pandas.container.DataFrameSchema)

#### *classmethod* empty(\*\_args) → [DataFrame](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.typing.DataFrame.html#pandera.typing.DataFrame)[Self]

Create an empty DataFrame with the schema of this model.

#### *classmethod* example(\*\*kwargs) → DataFrameBase[Self]

Generate an example of this data model specification.

#### *classmethod* from_json(source)

Load a schema from JSON.

* **Parameters:**
  **source** – str, Path, or file stream with JSON content.
* **Returns:**
  the backend-specific schema object.

#### *classmethod* from_yaml(yaml_schema)

Load a schema from YAML.

* **Parameters:**
  **yaml_schema** – str, Path, or file stream with YAML content.
* **Returns:**
  the backend-specific schema object.

#### *classmethod* get_metadata() → dict | None

Provide metadata for columns and schema level

#### *classmethod* pydantic_validate(schema_model: Any) → [DataFrameModel](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.api.dataframe.model.DataFrameModel.html#pandera.api.dataframe.model.DataFrameModel)

Verify that the input is a compatible dataframe model.

#### *classmethod* strategy(\*\*kwargs)

Create a data synthesis strategy.

#### *classmethod* to_json(target: PathLike | None = None, \*\*kwargs)

Convert this model’s schema to JSON.

#### *classmethod* to_json_schema()

Serialize schema metadata into json-schema format.

* **Parameters:**
  **dataframe_schema** – schema to write to json-schema format.

#### NOTE
This function is currently does not fully specify a pandera schema,
and is primarily used internally to render OpenAPI docs via the
FastAPI integration.

#### *classmethod* to_schema() → TSchema

Create `DataFrameSchema` from the `DataFrameModel`.

#### *classmethod* to_yaml(stream: PathLike | None = None)

Convert this model’s schema to YAML.

#### *classmethod* validate(check_obj: DataFrame, head: int | None = None, tail: int | None = None, sample: int | None = None, random_state: int | None = None, lazy: bool = False, inplace: bool = False) → [DataFrame](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.typing.DataFrame.html#pandera.typing.DataFrame)[Self]

Validate a DataFrame based on the schema specification.

* **Parameters:**
  * **check_obj** (*pd.DataFrame*) – the dataframe to be validated.
  * **head** – validate the first n rows. Rows overlapping with tail or
    sample are de-duplicated.
  * **tail** – validate the last n rows. Rows overlapping with head or
    sample are de-duplicated.
  * **sample** – validate a random sample of n rows. Rows overlapping
    with head or tail are de-duplicated.
  * **random_state** – random seed for the `sample` argument.
  * **lazy** – if True, lazily evaluates dataframe against all validation
    checks and raises a `SchemaErrors`. Otherwise, raise
    `SchemaError` as soon as one occurs.
  * **inplace** – if True, applies coercion to the object of validation,
    otherwise creates a copy of the data.
* **Returns:**
  validated `DataFrame`
* **Raises:**
  **SchemaError** – when `DataFrame` violates built-in or custom
  checks.

### *class* TwoLoadingsSchema(\*args, \*\*kwargs)

Bases: [`TwoComponentSchema`](#acore.types.exploratory_analysis.TwoComponentSchema)

Schema for the PCA loadings DataFrame.

#### value *: float | None* *= 'value'*

#### *class* Config

Bases: [`Config`](#acore.types.exploratory_analysis.TwoLoadingsSchema.Config)

#### add_missing_columns *: bool* *= False*

add columns to dataframe if they are missing

#### coerce *: bool* *= False*

coerce types of all schema components

#### description *: str | None* *= None*

arbitrary textual description

#### drop_invalid_rows *: bool* *= False*

drop invalid rows on validation

#### dtype *: str | type | [DataType](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.dtypes.DataType.html#pandera.dtypes.DataType) | ExtensionDtype | [dtype](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype) | None* *= None*

datatype of the dataframe. This overrides the data types specified in
any of the fields.

#### from_format *: Literal[Formats.csv] | Literal[Formats.dict] | Literal[Formats.json] | Literal[Formats.feather] | Literal[Formats.parquet] | Literal[Formats.pickle] | Literal[Formats.json_normalize] | Callable | None* *= None*

data format before validation. This option only applies to
schemas used in the context of the pandera type constructor
`pa.typing.DataFrame[Schema](data)`. If None, assumes a data structure
compatible with the `pandas.DataFrame` constructor.

#### from_format_kwargs *: dict[str, Any] | None* *= None*

a dictionary keyword arguments to pass into the reader function that
converts the object of type `from_format` to a pandera-validate-able
data structure. The reader function is implemented in the pandera.typing
generic types via the `from_format` and `to_format` methods.

#### metadata *: dict | None* *= None*

a dictionary object to store key-value data at schema level

#### multiindex_coerce *: bool* *= False*

coerce types of all MultiIndex components

#### multiindex_name *: str | None* *= None*

name of multiindex

#### multiindex_ordered *: bool* *= True*

validate MultiIndex in order

#### multiindex_strict *: bool | Literal['filter']* *= False*

make sure all specified columns are in validated MultiIndex -
if `"filter"`, removes indexes not specified in the schema

#### multiindex_unique *= None*

make sure the MultiIndex is unique along the list of columns

#### name *: str | None* *= 'TwoLoadingsSchema'*

name of schema

#### on_missing_columns *: str | None* *= None*

action to take when an optional column is missing from the dataframe.
`None` (default) passes silently, `"warn"` emits a
`SchemaWarning`. Overridable per field via
`Field(on_missing=...)`.

#### ordered *: bool* *= False*

validate columns order

#### strict *: bool | Literal['filter']* *= False*

make sure all specified columns are in the validated dataframe -
if `"filter"`, removes columns not specified in the schema

#### title *: str | None* *= None*

human-readable label for schema

#### to_format *: Literal[Formats.csv] | Literal[Formats.dict] | Literal[Formats.json] | Literal[Formats.feather] | Literal[Formats.parquet] | Literal[Formats.pickle] | Literal[Formats.json_normalize] | Callable | None* *= None*

data format to serialize into after validation. This option only applies
to  schemas used in the context of the pandera type constructor
`pa.typing.DataFrame[Schema](data)`. If None, returns a dataframe.

#### to_format_buffer *: str | Callable | None* *= None*

Buffer to be provided when to_format is a custom callable. See docs for
example of how to implement an example of a to format function.

#### to_format_kwargs *: dict[str, Any] | None* *= None*

a dictionary keyword arguments to pass into the writer function that
converts the pandera-validate-able object to type `to_format`.
The writer function is implemented in the pandera.typing
generic types via the `from_format` and `to_format` methods.

#### unique *: str | list[str] | None* *= None*

make sure certain column combinations are unique

#### unique_column_names *: bool* *= False*

make sure dataframe column names are unique

#### *classmethod* build_schema_(\*\*kwargs) → [DataFrameSchema](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.api.pandas.container.DataFrameSchema.html#pandera.api.pandas.container.DataFrameSchema)

#### *classmethod* empty(\*\_args) → [DataFrame](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.typing.DataFrame.html#pandera.typing.DataFrame)[Self]

Create an empty DataFrame with the schema of this model.

#### *classmethod* example(\*\*kwargs) → DataFrameBase[Self]

Generate an example of this data model specification.

#### *classmethod* from_json(source)

Load a schema from JSON.

* **Parameters:**
  **source** – str, Path, or file stream with JSON content.
* **Returns:**
  the backend-specific schema object.

#### *classmethod* from_yaml(yaml_schema)

Load a schema from YAML.

* **Parameters:**
  **yaml_schema** – str, Path, or file stream with YAML content.
* **Returns:**
  the backend-specific schema object.

#### *classmethod* get_metadata() → dict | None

Provide metadata for columns and schema level

#### group *: str | None* *= 'group'*

#### *classmethod* pydantic_validate(schema_model: Any) → [DataFrameModel](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.api.dataframe.model.DataFrameModel.html#pandera.api.dataframe.model.DataFrameModel)

Verify that the input is a compatible dataframe model.

#### *classmethod* strategy(\*\*kwargs)

Create a data synthesis strategy.

#### *classmethod* to_json(target: PathLike | None = None, \*\*kwargs)

Convert this model’s schema to JSON.

#### *classmethod* to_json_schema()

Serialize schema metadata into json-schema format.

* **Parameters:**
  **dataframe_schema** – schema to write to json-schema format.

#### NOTE
This function is currently does not fully specify a pandera schema,
and is primarily used internally to render OpenAPI docs via the
FastAPI integration.

#### *classmethod* to_schema() → TSchema

Create `DataFrameSchema` from the `DataFrameModel`.

#### *classmethod* to_yaml(stream: PathLike | None = None)

Convert this model’s schema to YAML.

#### *classmethod* validate(check_obj: DataFrame, head: int | None = None, tail: int | None = None, sample: int | None = None, random_state: int | None = None, lazy: bool = False, inplace: bool = False) → [DataFrame](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.typing.DataFrame.html#pandera.typing.DataFrame)[Self]

Validate a DataFrame based on the schema specification.

* **Parameters:**
  * **check_obj** (*pd.DataFrame*) – the dataframe to be validated.
  * **head** – validate the first n rows. Rows overlapping with tail or
    sample are de-duplicated.
  * **tail** – validate the last n rows. Rows overlapping with head or
    sample are de-duplicated.
  * **sample** – validate a random sample of n rows. Rows overlapping
    with head or tail are de-duplicated.
  * **random_state** – random seed for the `sample` argument.
  * **lazy** – if True, lazily evaluates dataframe against all validation
    checks and raises a `SchemaErrors`. Otherwise, raise
    `SchemaError` as soon as one occurs.
  * **inplace** – if True, applies coercion to the object of validation,
    otherwise creates a copy of the data.
* **Returns:**
  validated `DataFrame`
* **Raises:**
  **SchemaError** – when `DataFrame` violates built-in or custom
  checks.

#### x *: float* *= 'x'*

#### y *: float* *= 'y'*

### *class* AnnotationResult(, x_title: str, y_title: str, group: str | None = None)

Bases: [`BaseModel`](https://pydantic.dev/docs/validation/latest/api/pydantic/base_model/#pydantic.BaseModel)

Represents the annotation result from exploratory analysis.

#### x_title *: str*

#### y_title *: str*

#### group *: str | None*

#### model_config *= {}*

Configuration for the model, should be a dictionary conforming to [ConfigDict][pydantic.config.ConfigDict].

## acore.types.permutation_test module

### *class* PermutationResult(, metric: str | Callable | None = None, observed_statistic: Any = None, p_value: float)

Bases: [`BaseModel`](https://pydantic.dev/docs/validation/latest/api/pydantic/base_model/#pydantic.BaseModel)

Schema for the output of a permutation test.

#### metric *: str | Callable | None*

#### observed_statistic *: Any*

#### p_value *: float*

#### model_config *= {}*

Configuration for the model, should be a dictionary conforming to [ConfigDict][pydantic.config.ConfigDict].
