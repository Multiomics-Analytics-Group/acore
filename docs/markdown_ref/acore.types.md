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

#### group1 *: [str](https://docs.python.org/3/library/stdtypes.html#str)* *= 'group1'*

#### group2 *: [str](https://docs.python.org/3/library/stdtypes.html#str)* *= 'group2'*

#### mean_group1 *: [float](https://docs.python.org/3/library/functions.html#float)* *= 'mean(group1)'*

#### std_group1 *: [float](https://docs.python.org/3/library/functions.html#float)* *= 'std(group1)'*

#### mean_group2 *: [float](https://docs.python.org/3/library/functions.html#float)* *= 'mean(group2)'*

#### std_group2 *: [float](https://docs.python.org/3/library/functions.html#float)* *= 'std(group2)'*

#### t_statistics *: [float](https://docs.python.org/3/library/functions.html#float)* *= 'T-Statistics'*

#### pvalue *: [float](https://docs.python.org/3/library/functions.html#float)* *= 'pvalue'*

#### log2FC *: [float](https://docs.python.org/3/library/functions.html#float)* *= 'log2FC'*

#### FC *: [float](https://docs.python.org/3/library/functions.html#float)* *= 'FC'*

#### padj *: [float](https://docs.python.org/3/library/functions.html#float)* *= 'padj'*

#### correction *: [str](https://docs.python.org/3/library/stdtypes.html#str)* *= 'correction'*

#### rejected *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= 'rejected'*

#### neg_log10_p_value *: [float](https://docs.python.org/3/library/functions.html#float)* *= '-log10 pvalue'*

#### Method *: [str](https://docs.python.org/3/library/stdtypes.html#str)* *= 'Method'*

#### *class* Config

Bases: [`BaseConfig`](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.api.pandas.model_config.BaseConfig.html#pandera.api.pandas.model_config.BaseConfig)

#### add_missing_columns *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= False*

add columns to dataframe if they are missing

#### coerce *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= False*

coerce types of all schema components

#### description *: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

arbitrary textual description

#### drop_invalid_rows *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= False*

drop invalid rows on validation

#### dtype *: PandasDtypeInputTypes | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

datatype of the dataframe. This overrides the data types specified in
any of the fields.

#### from_format *: Format | Callable | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

data format before validation. This option only applies to
schemas used in the context of the pandera type constructor
`pa.typing.DataFrame[Schema](data)`. If None, assumes a data structure
compatible with the `pandas.DataFrame` constructor.

#### from_format_kwargs *: [dict](https://docs.python.org/3/library/stdtypes.html#dict)[[str](https://docs.python.org/3/library/stdtypes.html#str), Any] | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

a dictionary keyword arguments to pass into the reader function that
converts the object of type `from_format` to a pandera-validate-able
data structure. The reader function is implemented in the pandera.typing
generic types via the `from_format` and `to_format` methods.

#### metadata *: [dict](https://docs.python.org/3/library/stdtypes.html#dict) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

a dictionary object to store key-value data at schema level

#### multiindex_coerce *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= False*

coerce types of all MultiIndex components

#### multiindex_name *: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

name of multiindex

#### multiindex_ordered *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= True*

validate MultiIndex in order

#### multiindex_strict *: StrictType* *= False*

make sure all specified columns are in validated MultiIndex -
if `"filter"`, removes indexes not specified in the schema

#### multiindex_unique *= None*

make sure the MultiIndex is unique along the list of columns

#### name *: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None)* *= 'AnovaSchema'*

name of schema

#### ordered *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= False*

validate columns order

#### strict *: StrictType* *= False*

make sure all specified columns are in the validated dataframe -
if `"filter"`, removes columns not specified in the schema

#### title *: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

human-readable label for schema

#### to_format *: Format | Callable | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

data format to serialize into after validation. This option only applies
to  schemas used in the context of the pandera type constructor
`pa.typing.DataFrame[Schema](data)`. If None, returns a dataframe.

#### to_format_buffer *: [str](https://docs.python.org/3/library/stdtypes.html#str) | Callable | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

Buffer to be provided when to_format is a custom callable. See docs for
example of how to implement an example of a to format function.

#### to_format_kwargs *: [dict](https://docs.python.org/3/library/stdtypes.html#dict)[[str](https://docs.python.org/3/library/stdtypes.html#str), Any] | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

a dictionary keyword arguments to pass into the writer function that
converts the pandera-validate-able object to type `to_format`.
The writer function is implemented in the pandera.typing
generic types via the `from_format` and `to_format` methods.

#### unique *: [str](https://docs.python.org/3/library/stdtypes.html#str) | [list](https://docs.python.org/3/library/stdtypes.html#list)[[str](https://docs.python.org/3/library/stdtypes.html#str)] | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

make sure certain column combinations are unique

#### unique_column_names *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= False*

make sure dataframe column names are unique

#### *classmethod* build_schema_(\*\*kwargs) → [DataFrameSchema](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.api.pandas.container.DataFrameSchema.html#pandera.api.pandas.container.DataFrameSchema)

#### *classmethod* empty(\*\_args) → [DataFrame](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.typing.DataFrame.html#pandera.typing.DataFrame)[[Self](https://docs.python.org/3/library/typing.html#typing.Self)]

Create an empty DataFrame with the schema of this model.

#### *classmethod* example(\*\*kwargs) → DataFrameBase[TDataFrameModel]

Generate an example of a particular size.

* **Parameters:**
  **size** – number of elements in the generated DataFrame.
* **Returns:**
  DataFrame object.

#### *classmethod* get_metadata() → [dict](https://docs.python.org/3/library/stdtypes.html#dict) | [None](https://docs.python.org/3/library/constants.html#None)

Provide metadata for columns and schema level

#### *classmethod* pydantic_validate(schema_model: [Any](https://docs.python.org/3/library/typing.html#typing.Any)) → [DataFrameModel](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.api.dataframe.model.DataFrameModel.html#pandera.api.dataframe.model.DataFrameModel)

Verify that the input is a compatible dataframe model.

#### *classmethod* strategy(\*\*kwargs)

Create a `hypothesis` strategy for generating a DataFrame.

* **Parameters:**
  * **size** – number of elements to generate
  * **n_regex_columns** – number of regex columns to generate.
* **Returns:**
  a strategy that generates DataFrame objects.

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

#### *classmethod* to_yaml(stream: [PathLike](https://docs.python.org/3/library/os.html#os.PathLike) | [None](https://docs.python.org/3/library/constants.html#None) = None)

Convert Schema to yaml using io.to_yaml.

#### *classmethod* validate(check_obj: DataFrame, head: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None) = None, tail: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None) = None, sample: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None) = None, random_state: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None) = None, lazy: [bool](https://docs.python.org/3/library/functions.html#bool) = False, inplace: [bool](https://docs.python.org/3/library/functions.html#bool) = False) → [DataFrame](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.typing.DataFrame.html#pandera.typing.DataFrame)[[Self](https://docs.python.org/3/library/typing.html#typing.Self)]

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

#### t_statistics *: [float](https://docs.python.org/3/library/functions.html#float)* *= 'posthoc T-Statistics'*

#### posthoc_pvalue *: [float](https://docs.python.org/3/library/functions.html#float)* *= 'posthoc pvalue'*

#### f_statistics *: [float](https://docs.python.org/3/library/functions.html#float)* *= 'F-statistics'*

#### posthoc_padj *: [float](https://docs.python.org/3/library/functions.html#float)* *= 'posthoc padj'*

#### posthoc_paired *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= 'posthoc Paired'*

#### posthoc_parametric *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= 'posthoc Parametric'*

#### posthoc_dof *: [float](https://docs.python.org/3/library/functions.html#float)* *= 'posthoc dof'*

#### posthoc_tail *: [str](https://docs.python.org/3/library/stdtypes.html#str)* *= 'posthoc tail'*

#### posthoc_BF10 *: [str](https://docs.python.org/3/library/stdtypes.html#str)* *= 'posthoc BF10'*

#### posthoc_effsize *: [float](https://docs.python.org/3/library/functions.html#float)* *= 'posthoc effsize'*

#### efftype *: [str](https://docs.python.org/3/library/stdtypes.html#str)* *= 'efftype'*

#### *class* Config

Bases: [`Config`](#acore.types.exploratory_analysis.TwoLoadingsSchema.Config)

#### add_missing_columns *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= False*

add columns to dataframe if they are missing

#### coerce *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= False*

coerce types of all schema components

#### description *: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

arbitrary textual description

#### drop_invalid_rows *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= False*

drop invalid rows on validation

#### dtype *: PandasDtypeInputTypes | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

datatype of the dataframe. This overrides the data types specified in
any of the fields.

#### from_format *: Format | Callable | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

data format before validation. This option only applies to
schemas used in the context of the pandera type constructor
`pa.typing.DataFrame[Schema](data)`. If None, assumes a data structure
compatible with the `pandas.DataFrame` constructor.

#### from_format_kwargs *: [dict](https://docs.python.org/3/library/stdtypes.html#dict)[[str](https://docs.python.org/3/library/stdtypes.html#str), Any] | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

a dictionary keyword arguments to pass into the reader function that
converts the object of type `from_format` to a pandera-validate-able
data structure. The reader function is implemented in the pandera.typing
generic types via the `from_format` and `to_format` methods.

#### metadata *: [dict](https://docs.python.org/3/library/stdtypes.html#dict) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

a dictionary object to store key-value data at schema level

#### multiindex_coerce *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= False*

coerce types of all MultiIndex components

#### multiindex_name *: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

name of multiindex

#### multiindex_ordered *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= True*

validate MultiIndex in order

#### multiindex_strict *: StrictType* *= False*

make sure all specified columns are in validated MultiIndex -
if `"filter"`, removes indexes not specified in the schema

#### multiindex_unique *= None*

make sure the MultiIndex is unique along the list of columns

#### name *: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None)* *= 'AnovaSchemaMultiGroup'*

name of schema

#### ordered *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= False*

validate columns order

#### strict *: StrictType* *= False*

make sure all specified columns are in the validated dataframe -
if `"filter"`, removes columns not specified in the schema

#### title *: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

human-readable label for schema

#### to_format *: Format | Callable | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

data format to serialize into after validation. This option only applies
to  schemas used in the context of the pandera type constructor
`pa.typing.DataFrame[Schema](data)`. If None, returns a dataframe.

#### to_format_buffer *: [str](https://docs.python.org/3/library/stdtypes.html#str) | Callable | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

Buffer to be provided when to_format is a custom callable. See docs for
example of how to implement an example of a to format function.

#### to_format_kwargs *: [dict](https://docs.python.org/3/library/stdtypes.html#dict)[[str](https://docs.python.org/3/library/stdtypes.html#str), Any] | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

a dictionary keyword arguments to pass into the writer function that
converts the pandera-validate-able object to type `to_format`.
The writer function is implemented in the pandera.typing
generic types via the `from_format` and `to_format` methods.

#### unique *: [str](https://docs.python.org/3/library/stdtypes.html#str) | [list](https://docs.python.org/3/library/stdtypes.html#list)[[str](https://docs.python.org/3/library/stdtypes.html#str)] | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

make sure certain column combinations are unique

#### unique_column_names *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= False*

make sure dataframe column names are unique

#### FC *: [float](https://docs.python.org/3/library/functions.html#float)* *= 'FC'*

#### Method *: [str](https://docs.python.org/3/library/stdtypes.html#str)* *= 'Method'*

#### *classmethod* build_schema_(\*\*kwargs) → [DataFrameSchema](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.api.pandas.container.DataFrameSchema.html#pandera.api.pandas.container.DataFrameSchema)

#### correction *: [str](https://docs.python.org/3/library/stdtypes.html#str)* *= 'correction'*

#### *classmethod* empty(\*\_args) → [DataFrame](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.typing.DataFrame.html#pandera.typing.DataFrame)[[Self](https://docs.python.org/3/library/typing.html#typing.Self)]

Create an empty DataFrame with the schema of this model.

#### *classmethod* example(\*\*kwargs) → DataFrameBase[TDataFrameModel]

Generate an example of a particular size.

* **Parameters:**
  **size** – number of elements in the generated DataFrame.
* **Returns:**
  DataFrame object.

#### *classmethod* get_metadata() → [dict](https://docs.python.org/3/library/stdtypes.html#dict) | [None](https://docs.python.org/3/library/constants.html#None)

Provide metadata for columns and schema level

#### group1 *: [str](https://docs.python.org/3/library/stdtypes.html#str)* *= 'group1'*

#### group2 *: [str](https://docs.python.org/3/library/stdtypes.html#str)* *= 'group2'*

#### log2FC *: [float](https://docs.python.org/3/library/functions.html#float)* *= 'log2FC'*

#### mean_group1 *: [float](https://docs.python.org/3/library/functions.html#float)* *= 'mean(group1)'*

#### mean_group2 *: [float](https://docs.python.org/3/library/functions.html#float)* *= 'mean(group2)'*

#### neg_log10_p_value *: [float](https://docs.python.org/3/library/functions.html#float)* *= '-log10 pvalue'*

#### padj *: [float](https://docs.python.org/3/library/functions.html#float)* *= 'padj'*

#### pvalue *: [float](https://docs.python.org/3/library/functions.html#float)* *= 'pvalue'*

#### *classmethod* pydantic_validate(schema_model: [Any](https://docs.python.org/3/library/typing.html#typing.Any)) → [DataFrameModel](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.api.dataframe.model.DataFrameModel.html#pandera.api.dataframe.model.DataFrameModel)

Verify that the input is a compatible dataframe model.

#### rejected *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= 'rejected'*

#### std_group1 *: [float](https://docs.python.org/3/library/functions.html#float)* *= 'std(group1)'*

#### std_group2 *: [float](https://docs.python.org/3/library/functions.html#float)* *= 'std(group2)'*

#### *classmethod* strategy(\*\*kwargs)

Create a `hypothesis` strategy for generating a DataFrame.

* **Parameters:**
  * **size** – number of elements to generate
  * **n_regex_columns** – number of regex columns to generate.
* **Returns:**
  a strategy that generates DataFrame objects.

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

#### *classmethod* to_yaml(stream: [PathLike](https://docs.python.org/3/library/os.html#os.PathLike) | [None](https://docs.python.org/3/library/constants.html#None) = None)

Convert Schema to yaml using io.to_yaml.

#### *classmethod* validate(check_obj: DataFrame, head: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None) = None, tail: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None) = None, sample: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None) = None, random_state: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None) = None, lazy: [bool](https://docs.python.org/3/library/functions.html#bool) = False, inplace: [bool](https://docs.python.org/3/library/functions.html#bool) = False) → [DataFrame](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.typing.DataFrame.html#pandera.typing.DataFrame)[[Self](https://docs.python.org/3/library/typing.html#typing.Self)]

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

#### t_statistics *: [float](https://docs.python.org/3/library/functions.html#float)* *= 'posthoc T-Statistics'*

#### posthoc_pvalue *: [float](https://docs.python.org/3/library/functions.html#float)* *= 'posthoc pvalue'*

#### coef *: [float](https://docs.python.org/3/library/functions.html#float)* *= 'coef'*

#### std_err *: [float](https://docs.python.org/3/library/functions.html#float)* *= 'std err'*

#### conf_int_low *: [float](https://docs.python.org/3/library/functions.html#float)* *= 'Conf. Int. Low'*

#### conf_int_upp *: [float](https://docs.python.org/3/library/functions.html#float)* *= 'Conf. Int. Upp.'*

#### f_statistics *: [float](https://docs.python.org/3/library/functions.html#float)* *= 'F-statistics'*

#### posthoc_padj *: [float](https://docs.python.org/3/library/functions.html#float)* *= 'posthoc padj'*

#### *class* Config

Bases: [`BaseConfig`](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.api.pandas.model_config.BaseConfig.html#pandera.api.pandas.model_config.BaseConfig)

#### add_missing_columns *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= False*

add columns to dataframe if they are missing

#### coerce *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= False*

coerce types of all schema components

#### description *: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

arbitrary textual description

#### drop_invalid_rows *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= False*

drop invalid rows on validation

#### dtype *: PandasDtypeInputTypes | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

datatype of the dataframe. This overrides the data types specified in
any of the fields.

#### from_format *: Format | Callable | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

data format before validation. This option only applies to
schemas used in the context of the pandera type constructor
`pa.typing.DataFrame[Schema](data)`. If None, assumes a data structure
compatible with the `pandas.DataFrame` constructor.

#### from_format_kwargs *: [dict](https://docs.python.org/3/library/stdtypes.html#dict)[[str](https://docs.python.org/3/library/stdtypes.html#str), Any] | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

a dictionary keyword arguments to pass into the reader function that
converts the object of type `from_format` to a pandera-validate-able
data structure. The reader function is implemented in the pandera.typing
generic types via the `from_format` and `to_format` methods.

#### metadata *: [dict](https://docs.python.org/3/library/stdtypes.html#dict) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

a dictionary object to store key-value data at schema level

#### multiindex_coerce *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= False*

coerce types of all MultiIndex components

#### multiindex_name *: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

name of multiindex

#### multiindex_ordered *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= True*

validate MultiIndex in order

#### multiindex_strict *: StrictType* *= False*

make sure all specified columns are in validated MultiIndex -
if `"filter"`, removes indexes not specified in the schema

#### multiindex_unique *= None*

make sure the MultiIndex is unique along the list of columns

#### name *: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None)* *= 'AncovaSchema'*

name of schema

#### ordered *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= False*

validate columns order

#### strict *: StrictType* *= False*

make sure all specified columns are in the validated dataframe -
if `"filter"`, removes columns not specified in the schema

#### title *: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

human-readable label for schema

#### to_format *: Format | Callable | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

data format to serialize into after validation. This option only applies
to  schemas used in the context of the pandera type constructor
`pa.typing.DataFrame[Schema](data)`. If None, returns a dataframe.

#### to_format_buffer *: [str](https://docs.python.org/3/library/stdtypes.html#str) | Callable | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

Buffer to be provided when to_format is a custom callable. See docs for
example of how to implement an example of a to format function.

#### to_format_kwargs *: [dict](https://docs.python.org/3/library/stdtypes.html#dict)[[str](https://docs.python.org/3/library/stdtypes.html#str), Any] | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

a dictionary keyword arguments to pass into the writer function that
converts the pandera-validate-able object to type `to_format`.
The writer function is implemented in the pandera.typing
generic types via the `from_format` and `to_format` methods.

#### unique *: [str](https://docs.python.org/3/library/stdtypes.html#str) | [list](https://docs.python.org/3/library/stdtypes.html#list)[[str](https://docs.python.org/3/library/stdtypes.html#str)] | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

make sure certain column combinations are unique

#### unique_column_names *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= False*

make sure dataframe column names are unique

#### *classmethod* build_schema_(\*\*kwargs) → [DataFrameSchema](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.api.pandas.container.DataFrameSchema.html#pandera.api.pandas.container.DataFrameSchema)

#### *classmethod* empty(\*\_args) → [DataFrame](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.typing.DataFrame.html#pandera.typing.DataFrame)[[Self](https://docs.python.org/3/library/typing.html#typing.Self)]

Create an empty DataFrame with the schema of this model.

#### *classmethod* example(\*\*kwargs) → DataFrameBase[TDataFrameModel]

Generate an example of a particular size.

* **Parameters:**
  **size** – number of elements in the generated DataFrame.
* **Returns:**
  DataFrame object.

#### *classmethod* get_metadata() → [dict](https://docs.python.org/3/library/stdtypes.html#dict) | [None](https://docs.python.org/3/library/constants.html#None)

Provide metadata for columns and schema level

#### *classmethod* pydantic_validate(schema_model: [Any](https://docs.python.org/3/library/typing.html#typing.Any)) → [DataFrameModel](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.api.dataframe.model.DataFrameModel.html#pandera.api.dataframe.model.DataFrameModel)

Verify that the input is a compatible dataframe model.

#### *classmethod* strategy(\*\*kwargs)

Create a `hypothesis` strategy for generating a DataFrame.

* **Parameters:**
  * **size** – number of elements to generate
  * **n_regex_columns** – number of regex columns to generate.
* **Returns:**
  a strategy that generates DataFrame objects.

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

#### *classmethod* to_yaml(stream: [PathLike](https://docs.python.org/3/library/os.html#os.PathLike) | [None](https://docs.python.org/3/library/constants.html#None) = None)

Convert Schema to yaml using io.to_yaml.

#### *classmethod* validate(check_obj: DataFrame, head: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None) = None, tail: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None) = None, sample: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None) = None, random_state: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None) = None, lazy: [bool](https://docs.python.org/3/library/functions.html#bool) = False, inplace: [bool](https://docs.python.org/3/library/functions.html#bool) = False) → [DataFrame](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.typing.DataFrame.html#pandera.typing.DataFrame)[[Self](https://docs.python.org/3/library/typing.html#typing.Self)]

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

#### terms *: [str](https://docs.python.org/3/library/stdtypes.html#str)* *= 'terms'*

#### identifiers *: [str](https://docs.python.org/3/library/stdtypes.html#str)* *= 'identifiers'*

#### foreground *: [int](https://docs.python.org/3/library/functions.html#int)* *= 'foreground'*

#### background *: [int](https://docs.python.org/3/library/functions.html#int)* *= 'background'*

#### foreground_pop *: [int](https://docs.python.org/3/library/functions.html#int)* *= 'foreground_pop'*

#### background_pop *: [int](https://docs.python.org/3/library/functions.html#int)* *= 'background_pop'*

#### pvalue *: [float](https://docs.python.org/3/library/functions.html#float)* *= 'pvalue'*

#### padj *: [float](https://docs.python.org/3/library/functions.html#float)* *= 'padj'*

#### rejected *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= 'rejected'*

#### direction *: [str](https://docs.python.org/3/library/stdtypes.html#str)* *= 'direction'*

#### comparison *: [str](https://docs.python.org/3/library/stdtypes.html#str)* *= 'comparison'*

#### *class* Config

Bases: [`BaseConfig`](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.api.pandas.model_config.BaseConfig.html#pandera.api.pandas.model_config.BaseConfig)

#### add_missing_columns *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= False*

add columns to dataframe if they are missing

#### coerce *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= False*

coerce types of all schema components

#### description *: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

arbitrary textual description

#### drop_invalid_rows *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= False*

drop invalid rows on validation

#### dtype *: PandasDtypeInputTypes | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

datatype of the dataframe. This overrides the data types specified in
any of the fields.

#### from_format *: Format | Callable | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

data format before validation. This option only applies to
schemas used in the context of the pandera type constructor
`pa.typing.DataFrame[Schema](data)`. If None, assumes a data structure
compatible with the `pandas.DataFrame` constructor.

#### from_format_kwargs *: [dict](https://docs.python.org/3/library/stdtypes.html#dict)[[str](https://docs.python.org/3/library/stdtypes.html#str), Any] | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

a dictionary keyword arguments to pass into the reader function that
converts the object of type `from_format` to a pandera-validate-able
data structure. The reader function is implemented in the pandera.typing
generic types via the `from_format` and `to_format` methods.

#### metadata *: [dict](https://docs.python.org/3/library/stdtypes.html#dict) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

a dictionary object to store key-value data at schema level

#### multiindex_coerce *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= False*

coerce types of all MultiIndex components

#### multiindex_name *: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

name of multiindex

#### multiindex_ordered *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= True*

validate MultiIndex in order

#### multiindex_strict *: StrictType* *= False*

make sure all specified columns are in validated MultiIndex -
if `"filter"`, removes indexes not specified in the schema

#### multiindex_unique *= None*

make sure the MultiIndex is unique along the list of columns

#### name *: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None)* *= 'EnrichmentAnalysisSchema'*

name of schema

#### ordered *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= False*

validate columns order

#### strict *: StrictType* *= False*

make sure all specified columns are in the validated dataframe -
if `"filter"`, removes columns not specified in the schema

#### title *: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

human-readable label for schema

#### to_format *: Format | Callable | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

data format to serialize into after validation. This option only applies
to  schemas used in the context of the pandera type constructor
`pa.typing.DataFrame[Schema](data)`. If None, returns a dataframe.

#### to_format_buffer *: [str](https://docs.python.org/3/library/stdtypes.html#str) | Callable | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

Buffer to be provided when to_format is a custom callable. See docs for
example of how to implement an example of a to format function.

#### to_format_kwargs *: [dict](https://docs.python.org/3/library/stdtypes.html#dict)[[str](https://docs.python.org/3/library/stdtypes.html#str), Any] | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

a dictionary keyword arguments to pass into the writer function that
converts the pandera-validate-able object to type `to_format`.
The writer function is implemented in the pandera.typing
generic types via the `from_format` and `to_format` methods.

#### unique *: [str](https://docs.python.org/3/library/stdtypes.html#str) | [list](https://docs.python.org/3/library/stdtypes.html#list)[[str](https://docs.python.org/3/library/stdtypes.html#str)] | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

make sure certain column combinations are unique

#### unique_column_names *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= False*

make sure dataframe column names are unique

#### *classmethod* build_schema_(\*\*kwargs) → [DataFrameSchema](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.api.pandas.container.DataFrameSchema.html#pandera.api.pandas.container.DataFrameSchema)

#### *classmethod* empty(\*\_args) → [DataFrame](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.typing.DataFrame.html#pandera.typing.DataFrame)[[Self](https://docs.python.org/3/library/typing.html#typing.Self)]

Create an empty DataFrame with the schema of this model.

#### *classmethod* example(\*\*kwargs) → DataFrameBase[TDataFrameModel]

Generate an example of a particular size.

* **Parameters:**
  **size** – number of elements in the generated DataFrame.
* **Returns:**
  DataFrame object.

#### *classmethod* get_metadata() → [dict](https://docs.python.org/3/library/stdtypes.html#dict) | [None](https://docs.python.org/3/library/constants.html#None)

Provide metadata for columns and schema level

#### *classmethod* pydantic_validate(schema_model: [Any](https://docs.python.org/3/library/typing.html#typing.Any)) → [DataFrameModel](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.api.dataframe.model.DataFrameModel.html#pandera.api.dataframe.model.DataFrameModel)

Verify that the input is a compatible dataframe model.

#### *classmethod* strategy(\*\*kwargs)

Create a `hypothesis` strategy for generating a DataFrame.

* **Parameters:**
  * **size** – number of elements to generate
  * **n_regex_columns** – number of regex columns to generate.
* **Returns:**
  a strategy that generates DataFrame objects.

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

#### *classmethod* to_yaml(stream: [PathLike](https://docs.python.org/3/library/os.html#os.PathLike) | [None](https://docs.python.org/3/library/constants.html#None) = None)

Convert Schema to yaml using io.to_yaml.

#### *classmethod* validate(check_obj: DataFrame, head: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None) = None, tail: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None) = None, sample: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None) = None, random_state: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None) = None, lazy: [bool](https://docs.python.org/3/library/functions.html#bool) = False, inplace: [bool](https://docs.python.org/3/library/functions.html#bool) = False) → [DataFrame](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.typing.DataFrame.html#pandera.typing.DataFrame)[[Self](https://docs.python.org/3/library/typing.html#typing.Self)]

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

#### group *: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None)* *= 'group'*

#### x *: [float](https://docs.python.org/3/library/functions.html#float)* *= 'x'*

#### y *: [float](https://docs.python.org/3/library/functions.html#float)* *= 'y'*

#### *class* Config

Bases: [`BaseConfig`](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.api.pandas.model_config.BaseConfig.html#pandera.api.pandas.model_config.BaseConfig)

#### add_missing_columns *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= False*

add columns to dataframe if they are missing

#### coerce *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= False*

coerce types of all schema components

#### description *: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

arbitrary textual description

#### drop_invalid_rows *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= False*

drop invalid rows on validation

#### dtype *: PandasDtypeInputTypes | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

datatype of the dataframe. This overrides the data types specified in
any of the fields.

#### from_format *: Format | Callable | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

data format before validation. This option only applies to
schemas used in the context of the pandera type constructor
`pa.typing.DataFrame[Schema](data)`. If None, assumes a data structure
compatible with the `pandas.DataFrame` constructor.

#### from_format_kwargs *: [dict](https://docs.python.org/3/library/stdtypes.html#dict)[[str](https://docs.python.org/3/library/stdtypes.html#str), Any] | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

a dictionary keyword arguments to pass into the reader function that
converts the object of type `from_format` to a pandera-validate-able
data structure. The reader function is implemented in the pandera.typing
generic types via the `from_format` and `to_format` methods.

#### metadata *: [dict](https://docs.python.org/3/library/stdtypes.html#dict) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

a dictionary object to store key-value data at schema level

#### multiindex_coerce *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= False*

coerce types of all MultiIndex components

#### multiindex_name *: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

name of multiindex

#### multiindex_ordered *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= True*

validate MultiIndex in order

#### multiindex_strict *: StrictType* *= False*

make sure all specified columns are in validated MultiIndex -
if `"filter"`, removes indexes not specified in the schema

#### multiindex_unique *= None*

make sure the MultiIndex is unique along the list of columns

#### name *: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None)* *= 'TwoComponentSchema'*

name of schema

#### ordered *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= False*

validate columns order

#### strict *: StrictType* *= False*

make sure all specified columns are in the validated dataframe -
if `"filter"`, removes columns not specified in the schema

#### title *: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

human-readable label for schema

#### to_format *: Format | Callable | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

data format to serialize into after validation. This option only applies
to  schemas used in the context of the pandera type constructor
`pa.typing.DataFrame[Schema](data)`. If None, returns a dataframe.

#### to_format_buffer *: [str](https://docs.python.org/3/library/stdtypes.html#str) | Callable | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

Buffer to be provided when to_format is a custom callable. See docs for
example of how to implement an example of a to format function.

#### to_format_kwargs *: [dict](https://docs.python.org/3/library/stdtypes.html#dict)[[str](https://docs.python.org/3/library/stdtypes.html#str), Any] | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

a dictionary keyword arguments to pass into the writer function that
converts the pandera-validate-able object to type `to_format`.
The writer function is implemented in the pandera.typing
generic types via the `from_format` and `to_format` methods.

#### unique *: [str](https://docs.python.org/3/library/stdtypes.html#str) | [list](https://docs.python.org/3/library/stdtypes.html#list)[[str](https://docs.python.org/3/library/stdtypes.html#str)] | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

make sure certain column combinations are unique

#### unique_column_names *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= False*

make sure dataframe column names are unique

#### *classmethod* build_schema_(\*\*kwargs) → [DataFrameSchema](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.api.pandas.container.DataFrameSchema.html#pandera.api.pandas.container.DataFrameSchema)

#### *classmethod* empty(\*\_args) → [DataFrame](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.typing.DataFrame.html#pandera.typing.DataFrame)[[Self](https://docs.python.org/3/library/typing.html#typing.Self)]

Create an empty DataFrame with the schema of this model.

#### *classmethod* example(\*\*kwargs) → DataFrameBase[TDataFrameModel]

Generate an example of a particular size.

* **Parameters:**
  **size** – number of elements in the generated DataFrame.
* **Returns:**
  DataFrame object.

#### *classmethod* get_metadata() → [dict](https://docs.python.org/3/library/stdtypes.html#dict) | [None](https://docs.python.org/3/library/constants.html#None)

Provide metadata for columns and schema level

#### *classmethod* pydantic_validate(schema_model: [Any](https://docs.python.org/3/library/typing.html#typing.Any)) → [DataFrameModel](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.api.dataframe.model.DataFrameModel.html#pandera.api.dataframe.model.DataFrameModel)

Verify that the input is a compatible dataframe model.

#### *classmethod* strategy(\*\*kwargs)

Create a `hypothesis` strategy for generating a DataFrame.

* **Parameters:**
  * **size** – number of elements to generate
  * **n_regex_columns** – number of regex columns to generate.
* **Returns:**
  a strategy that generates DataFrame objects.

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

#### *classmethod* to_yaml(stream: [PathLike](https://docs.python.org/3/library/os.html#os.PathLike) | [None](https://docs.python.org/3/library/constants.html#None) = None)

Convert Schema to yaml using io.to_yaml.

#### *classmethod* validate(check_obj: DataFrame, head: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None) = None, tail: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None) = None, sample: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None) = None, random_state: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None) = None, lazy: [bool](https://docs.python.org/3/library/functions.html#bool) = False, inplace: [bool](https://docs.python.org/3/library/functions.html#bool) = False) → [DataFrame](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.typing.DataFrame.html#pandera.typing.DataFrame)[[Self](https://docs.python.org/3/library/typing.html#typing.Self)]

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

#### value *: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None)* *= 'value'*

#### *class* Config

Bases: [`Config`](#acore.types.exploratory_analysis.TwoLoadingsSchema.Config)

#### add_missing_columns *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= False*

add columns to dataframe if they are missing

#### coerce *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= False*

coerce types of all schema components

#### description *: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

arbitrary textual description

#### drop_invalid_rows *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= False*

drop invalid rows on validation

#### dtype *: PandasDtypeInputTypes | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

datatype of the dataframe. This overrides the data types specified in
any of the fields.

#### from_format *: Format | Callable | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

data format before validation. This option only applies to
schemas used in the context of the pandera type constructor
`pa.typing.DataFrame[Schema](data)`. If None, assumes a data structure
compatible with the `pandas.DataFrame` constructor.

#### from_format_kwargs *: [dict](https://docs.python.org/3/library/stdtypes.html#dict)[[str](https://docs.python.org/3/library/stdtypes.html#str), Any] | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

a dictionary keyword arguments to pass into the reader function that
converts the object of type `from_format` to a pandera-validate-able
data structure. The reader function is implemented in the pandera.typing
generic types via the `from_format` and `to_format` methods.

#### metadata *: [dict](https://docs.python.org/3/library/stdtypes.html#dict) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

a dictionary object to store key-value data at schema level

#### multiindex_coerce *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= False*

coerce types of all MultiIndex components

#### multiindex_name *: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

name of multiindex

#### multiindex_ordered *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= True*

validate MultiIndex in order

#### multiindex_strict *: StrictType* *= False*

make sure all specified columns are in validated MultiIndex -
if `"filter"`, removes indexes not specified in the schema

#### multiindex_unique *= None*

make sure the MultiIndex is unique along the list of columns

#### name *: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None)* *= 'TwoLoadingsSchema'*

name of schema

#### ordered *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= False*

validate columns order

#### strict *: StrictType* *= False*

make sure all specified columns are in the validated dataframe -
if `"filter"`, removes columns not specified in the schema

#### title *: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

human-readable label for schema

#### to_format *: Format | Callable | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

data format to serialize into after validation. This option only applies
to  schemas used in the context of the pandera type constructor
`pa.typing.DataFrame[Schema](data)`. If None, returns a dataframe.

#### to_format_buffer *: [str](https://docs.python.org/3/library/stdtypes.html#str) | Callable | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

Buffer to be provided when to_format is a custom callable. See docs for
example of how to implement an example of a to format function.

#### to_format_kwargs *: [dict](https://docs.python.org/3/library/stdtypes.html#dict)[[str](https://docs.python.org/3/library/stdtypes.html#str), Any] | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

a dictionary keyword arguments to pass into the writer function that
converts the pandera-validate-able object to type `to_format`.
The writer function is implemented in the pandera.typing
generic types via the `from_format` and `to_format` methods.

#### unique *: [str](https://docs.python.org/3/library/stdtypes.html#str) | [list](https://docs.python.org/3/library/stdtypes.html#list)[[str](https://docs.python.org/3/library/stdtypes.html#str)] | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

make sure certain column combinations are unique

#### unique_column_names *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= False*

make sure dataframe column names are unique

#### *classmethod* build_schema_(\*\*kwargs) → [DataFrameSchema](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.api.pandas.container.DataFrameSchema.html#pandera.api.pandas.container.DataFrameSchema)

#### *classmethod* empty(\*\_args) → [DataFrame](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.typing.DataFrame.html#pandera.typing.DataFrame)[[Self](https://docs.python.org/3/library/typing.html#typing.Self)]

Create an empty DataFrame with the schema of this model.

#### *classmethod* example(\*\*kwargs) → DataFrameBase[TDataFrameModel]

Generate an example of a particular size.

* **Parameters:**
  **size** – number of elements in the generated DataFrame.
* **Returns:**
  DataFrame object.

#### *classmethod* get_metadata() → [dict](https://docs.python.org/3/library/stdtypes.html#dict) | [None](https://docs.python.org/3/library/constants.html#None)

Provide metadata for columns and schema level

#### group *: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None)* *= 'group'*

#### *classmethod* pydantic_validate(schema_model: [Any](https://docs.python.org/3/library/typing.html#typing.Any)) → [DataFrameModel](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.api.dataframe.model.DataFrameModel.html#pandera.api.dataframe.model.DataFrameModel)

Verify that the input is a compatible dataframe model.

#### *classmethod* strategy(\*\*kwargs)

Create a `hypothesis` strategy for generating a DataFrame.

* **Parameters:**
  * **size** – number of elements to generate
  * **n_regex_columns** – number of regex columns to generate.
* **Returns:**
  a strategy that generates DataFrame objects.

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

#### *classmethod* to_yaml(stream: [PathLike](https://docs.python.org/3/library/os.html#os.PathLike) | [None](https://docs.python.org/3/library/constants.html#None) = None)

Convert Schema to yaml using io.to_yaml.

#### *classmethod* validate(check_obj: DataFrame, head: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None) = None, tail: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None) = None, sample: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None) = None, random_state: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None) = None, lazy: [bool](https://docs.python.org/3/library/functions.html#bool) = False, inplace: [bool](https://docs.python.org/3/library/functions.html#bool) = False) → [DataFrame](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.typing.DataFrame.html#pandera.typing.DataFrame)[[Self](https://docs.python.org/3/library/typing.html#typing.Self)]

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

#### x *: [float](https://docs.python.org/3/library/functions.html#float)* *= 'x'*

#### y *: [float](https://docs.python.org/3/library/functions.html#float)* *= 'y'*

### *class* AnnotationResult(\*, x_title: [str](https://docs.python.org/3/library/stdtypes.html#str), y_title: [str](https://docs.python.org/3/library/stdtypes.html#str), group: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None) = None)

Bases: [`BaseModel`](https://docs.pydantic.dev/latest/api/base_model/#pydantic.BaseModel)

Represents the annotation result from exploratory analysis.

#### x_title *: [str](https://docs.python.org/3/library/stdtypes.html#str)*

#### y_title *: [str](https://docs.python.org/3/library/stdtypes.html#str)*

#### group *: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None)*

#### model_config *: ClassVar[ConfigDict]* *= {}*

Configuration for the model, should be a dictionary conforming to [ConfigDict][pydantic.config.ConfigDict].
