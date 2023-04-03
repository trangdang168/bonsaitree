from druid_amish import *
import ast

def collect_node_dict(file_name, druid_df):
    logging.info("--------Collect node dicts!--------")
    po, sources, married_in, genotyped, mcras_df = amish_set_up()
    node_dict_df = druid_df[["real", "ca1", "ca2"]]
    node_dict_df["real"] = node_dict_df["real"].map(ast.literal_eval)
    # TODO FIX
    node_dict_df[] = node_dict_df.apply(lambda row: row["real"].update(get_node_dict_for_root(row["ca1"], po)) ,axis=0)
    node_dict_df.apply(lambda row: row["real"].update(get_node_dict_for_root(row["ca2"], po)) ,axis=0)
    
    node_dict_df.to_csv(file_name)

def main():
    druid_df = pd.read_csv(sys.argv[2])
    collect_node_dict(sys.argv[1], druid_df)

if __name__ == "__main__":
    main()